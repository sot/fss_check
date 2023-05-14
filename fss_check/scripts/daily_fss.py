# Licensed under a 3-clause BSD style license - see LICENSE.rst

import argparse
from datetime import timezone
from pathlib import Path

import astropy.units as u
import jinja2
import matplotlib
import numpy as np
from acdc.common import send_mail
from astropy.table import vstack
from cheta import fetch_eng as fetch
from cxotime import CxoTime, CxoTimeLike
from ska_helpers.logging import basic_logger
import yaml

import fss_check
from fss_check.check_fss import (
    get_large_pitch_roll_error_intervals,
    plot_delta_vs_pitch_roll,
    plot_pitch_for_data_with_large_errors,
    plot_pitch_roll_spm_mp_constraints,
    plot_roll_pitch_vs_time,
)
from fss_check.fss_utils import (
    add_pitch_roll_columns,
    get_fss_prim_data,
    get_fss_prim_data_cached,
)

logger = basic_logger("fss_check")

# Global config for convenience. This gets updated in main() after reading config file.
CONFIG = {}


def get_parser():
    parser = argparse.ArgumentParser(
        description="Daily Fine Sun Sensor monitoring and trending"
    )
    parser.add_argument("--out", type=str, default=".", help="Output data directory")
    parser.add_argument(
        "--config-dir",
        type=str,
        help="Config directory (default=current directory)",
    )
    parser.add_argument(
        "--days-recent",
        default=30,
        type=int,
        help="Number of days before --stop for recent plots",
    )
    parser.add_argument(
        "--days-table",
        default=90,
        type=int,
        help="Number of days before --stop for table of large pitch/roll errors",
    )
    parser.add_argument(
        "--days-long-term",
        default=360,
        type=int,
        help="Number of days before --stop for long term plots",
    )
    parser.add_argument(
        "--highlight-recent-days",
        type=float,
        default=30.0,
        help="Number of days to highlight in plots and table (days, default=30)",
    )
    parser.add_argument("--stop", type=str, help="Stop date (default=now)")
    parser.add_argument(
        "--use-maude",
        action="store_true",
        default=False,
        help="Use MAUDE to get most recent data",
    )
    parser.add_argument(
        "--email",
        action="append",
        dest="emails",
        default=[],
        help="Email address for notification (multiple allowed)",
    )
    parser.add_argument(
        "--level",
        default="INFO",
        help="Logging level (default=INFO)",
    )
    # Add --cache-data argument to cache data locally
    parser.add_argument(
        "--cache-data",
        action="store_true",
        default=False,
        help="Cache data locally in ./cache (mostly for testing)",
    )
    return parser


def main(args=None):
    data_dir = Path(fss_check.__file__).parent / "data"

    matplotlib.use("Agg")

    args = get_parser().parse_args(args)

    # Update global config with config file
    CONFIG.update(get_config(args.config_dir))

    if args.days_recent > args.days_long_term:
        raise ValueError("recent days must be <= long term days")
    if args.use_maude:
        fetch.data_source.set("cxc", "maude allow_subset=False")
    logger.setLevel(args.level.upper())

    stop = CxoTime(args.stop)
    starts = {
        "long_term": stop - args.days_long_term * u.day,
        "recent": stop - args.days_recent * u.day,
        "table": stop - args.days_table * u.day,
    }

    # Get data
    dats = get_datasets(args.cache_data, stop, starts)

    # Last available date in data
    date_last = CxoTime(dats["recent"]["times"][-1])
    date_last_local_fmt = get_formatted_local_time(date_last)

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # Classic plot of points with large > 2 deg pitch or roll errors near FSS FOV edge
    large_pr_err_sun, large_pr_err_no_sun = process_large_pitch_roll_errors(
        args, stop, starts, dats, outdir
    )

    pitch_roll_time_outfiles = process_pitch_roll_time_plots(
        dats["recent"], outdir, large_pr_err_sun, large_pr_err_no_sun
    )

    # 2 x 2 grid of pitch/roll error vs pitch/roll for data within SPM pitch limit
    plot_delta_vs_pitch_roll(
        dats["recent"],
        max_pitch=CONFIG["spm_pitch_limit"],
        err_lim=CONFIG["plot_delta_vs_pitch_roll"]["err_lim"],
        outfile=outdir / "delta_pitch_roll_vs_pitch_roll_recent.png",
    )

    plot_pitch_roll_spm_mp_constraints(
        dats["recent"],
        pitch_max=CONFIG["spm_pitch_limit"],
        err_caution=CONFIG["plot_pitch_roll_spm_mp_constraints"]["err_caution"],
        err_warning=CONFIG["plot_pitch_roll_spm_mp_constraints"]["err_warning"],
        outfile=outdir / "pitch_roll_spm_mp_constraints_recent.png",
    )

    # Write out the HTML report
    template = jinja2.Template((data_dir / "index.html").read_text())
    context = {
        "date_last": date_last.date[:-4],
        "date_last_local": date_last_local_fmt,
        "stop": stop.date[:-4],
        "days_long_term": args.days_long_term,
        "days_recent": args.days_recent,
        "days_table": args.days_table,
        "large_pitch_roll_error": large_pr_err_sun,
        "large_prerr_no_sun": large_pr_err_no_sun,
        "config": CONFIG,
        "axes": ["pitch", "roll"],
        "pitch_roll_time_outfiles": pitch_roll_time_outfiles,
    }
    txt = template.render(**context)
    (outdir / "index.html").write_text(txt)


def get_config(config_dir=None):
    """Load fss_check_config.yaml from config_dir (if defined) or package"""
    config_dir = Path(config_dir) if config_dir else Path(fss_check.__file__).parent
    path = config_dir / "fss_check_config.yml"
    logger.info(f"Loading config from {path}")
    with open(path) as fh:
        config = yaml.safe_load(fh)

    config["config_path"] = str(path.absolute())
    config["config_text"] = path.read_text()

    return config


def process_pitch_roll_time_plots(dat, outdir, large_pr_err_sun, large_pr_err_no_sun):
    logger.info("Making pitch roll detail plots")
    pitch_roll_time_outfiles = []
    pitch_roll_time_table = vstack(
        [
            large_pr_err_sun[: CONFIG["max_pitch_roll_time_plots"]],
            large_pr_err_no_sun[: CONFIG["max_pitch_roll_time_plots"]],
        ]
    )

    for interval in pitch_roll_time_table:
        if not interval["recent"]:
            continue
        start = CxoTime(interval["datestart"]) - 5 * u.min
        stop = CxoTime(interval["datestop"]) + 5 * u.min
        outfile = outdir / f"pitch_roll_vs_time_{interval['datestart']}.png"
        plot_roll_pitch_vs_time(
            dat,
            start,
            stop,
            pitch_max=CONFIG["spm_pitch_warning"],
            plot_errs=CONFIG["plot_roll_pitch_vs_time"]["plot_errs"],
            outfile=outfile,
            suptitle=interval["datestart"],
        )
        pitch_roll_time_outfiles.append(outfile.name)

    return pitch_roll_time_outfiles


def process_large_pitch_roll_errors(args, stop, starts, dats, outdir):
    logger.info("Finding large pitch/roll errors")
    for i_epoch, epoch in enumerate(starts):
        start = starts[epoch]
        dat = dats[epoch]
        for axis in ["pitch", "roll"]:
            show_legend = i_epoch == 0 and axis == "pitch"
            plot_pitch_for_data_with_large_errors(
                dat=dat,
                start=start,
                stop=stop,
                pitch_warning=CONFIG["spm_pitch_warning"],
                pitch_limit=CONFIG["spm_pitch_limit"],
                plot_pitch_min=CONFIG["plot_pitch_for_data_with_large_errors"][
                    "plot_pitch_min"
                ],
                outfile=outdir / f"pitch_bad_{axis}_{epoch}.png",
                axis=axis,
                show_legend=show_legend,
            )

    # Table of intervals of large (> 2 deg) pitch or roll errors
    large_pr_error_sun = get_large_pitch_roll_error_intervals(
        dats["table"],
        pitch_max=CONFIG["spm_pitch_warning"],
        err_min=CONFIG["get_large_pitch_roll_error_intervals"]["err_min"],
        dt_join=CONFIG["get_large_pitch_roll_error_intervals"]["dt_join"],
        sun_presence=True,
    )
    dt = stop - CxoTime(large_pr_error_sun["datestart"])
    large_pr_error_sun["recent"] = dt < args.highlight_recent_days * u.day

    # Table of intervals of large (> 2 deg) pitch or roll errors
    large_prerr_no_sun = get_large_pitch_roll_error_intervals(
        dats["table"],
        pitch_max=CONFIG["spm_pitch_limit"],
        err_min=CONFIG["get_large_pitch_roll_error_intervals"]["err_min"],
        dt_join=CONFIG["get_large_pitch_roll_error_intervals"]["dt_join"],
        sun_presence=False,
    )
    dt = stop - CxoTime(large_prerr_no_sun["datestart"])
    large_prerr_no_sun["recent"] = dt < args.highlight_recent_days * u.day

    return large_pr_error_sun, large_prerr_no_sun


def get_datasets(cache_data, stop, starts):
    logger.info(f"Processing {starts['long_term'].date} to {stop.date}")
    get_data_func = get_fss_prim_data_cached if cache_data else get_fss_prim_data
    dat = get_data_func(
        starts["long_term"],
        stop,
        exclude_intervals=CONFIG["exclude_intervals"],
    )
    dat = add_pitch_roll_columns(dat)

    # Extract a recent data subset
    dats = {"long_term": dat}
    for subset in ["recent", "table"]:
        i0 = np.searchsorted(dat["times"], starts[subset].cxcsec)
        dats[subset] = dat[i0:]
    return dats


def send_alert_email(opt, logger):
    subject = f"fss_check3: large pitch/roll error"
    lines = ["Long drop interval(s) found for the following perigee events:"]
    text = "\n".join(lines)
    send_mail(logger, opt, subject, text, __file__)


def get_formatted_local_time(date: CxoTimeLike):
    date = CxoTime(date)
    date_last_local = date.datetime.replace(tzinfo=timezone.utc).astimezone(tz=None)
    date_last_local_fmt = date_last_local.strftime("%Y %a %b %d %I:%M:%S %p %Z")
    return date_last_local_fmt


if __name__ == "__main__":
    main()
