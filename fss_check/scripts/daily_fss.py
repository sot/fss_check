# Licensed under a 3-clause BSD style license - see LICENSE.rst

import argparse
import shutil
from datetime import timezone
from pathlib import Path

import astropy.units as u
import jinja2
import matplotlib
import numpy as np
from cheta import fetch_eng as fetch
from cxotime import CxoTime, CxoTimeLike
from ska_helpers.logging import basic_logger

import fss_check
import fss_check.config
from fss_check.check_fss import (
    get_large_pitch_roll_error_intervals,
    plot_delta_vs_pitch_roll,
    plot_pitch_for_data_with_large_errors,
    plot_roll_pitch_vs_time,
)
from fss_check.config import CONFIG
from fss_check.fss_utils import (
    add_pitch_roll_columns,
    get_fss_prim_data,
    get_fss_prim_data_cached,
)

logger = basic_logger("fss_check")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Daily Fine Sun Sensor monitoring and trending"
    )
    parser.add_argument("--out", type=str, default=".", help="Output data directory")
    parser.add_argument(
        "--config-dir",
        type=str,
        default='.',
        help="Config directory (default=current directory)",
    )
    parser.add_argument(
        "--days-recent",
        default=10,
        type=int,
        help="Number of days before --stop for recent plots",
    )
    parser.add_argument(
        "--days-long-term",
        default=20,
        type=int,
        help="Number of days before --stop for long term plots",
    )
    parser.add_argument("--stop", type=str, help="Stop date (default=now)")
    parser.add_argument(
        "--use-maude",
        action="store_true",
        default=False,
        help="Use MAUDE to get most recent data",
    )
    # Add --level argument to set logging level
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
    parser.add_argument(
        "--highlight-recent-days",
        type=float,
        default=30.0,
        help="Number of days to highlight in plots and table (days, default=30)",
    )
    return parser


def main(args=None):
    data_dir = Path(fss_check.__file__).parent / "data"

    matplotlib.use("Agg")

    args = get_parser().parse_args(args)

    # TODO: improve this
    fss_check.config.CONFIG_DIR = args.config_dir

    if args.days_recent > args.days_long_term:
        raise ValueError("recent days must be <= long term days")
    if args.use_maude:
        fetch.data_source.set("cxc", "maude allow_subset=False")
    logger.setLevel(args.level.upper())

    stop = CxoTime(args.stop)
    starts = {
        "long_term": stop - args.days_long_term * u.day,
        "recent": stop - args.days_recent * u.day,
    }

    # Get data
    logger.info(f"Processing {starts['long_term'].date} to {stop.date}")
    get_data_func = get_fss_prim_data_cached if args.cache_data else get_fss_prim_data
    dat = get_data_func(starts["long_term"], stop)
    dat = add_pitch_roll_columns(dat)

    # Extract a recent data subset
    i0 = np.searchsorted(dat["times"], starts["recent"].cxcsec)
    dats = {"long_term": dat, "recent": dat[i0:]}

    # Last available date in data
    date_last = CxoTime(dat["times"][-1])
    date_last_local_fmt = get_formatted_local_time(date_last)

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # Classic plot of points with large > 2 deg pitch or roll errors
    for epoch in starts:
        start = starts[epoch]
        dat = dats[epoch]
        for axis in ["pitch", "roll"]:
            plot_pitch_for_data_with_large_errors(
                dat=dat,
                start=start,
                stop=stop,
                outfile=outdir / f"pitch_bad_{axis}_{epoch}.png",
                axis=axis,
            )

    # 2 x 2 grid of pitch/roll error vs pitch/roll
    plot_delta_vs_pitch_roll(
        dats["recent"], outfile=outdir / "delta_pitch_roll_vs_pitch_roll_recent.png"
    )

    # Table of intervals of large (> 2 deg) pitch or roll errors
    large_pitch_roll_error = get_large_pitch_roll_error_intervals(
        dats["recent"],
        pitch_max=CONFIG["spm_pitch_warning"],
        err_min=CONFIG["get_large_pitch_roll_error_intervals"]["err_min"],
        dt_join=CONFIG["get_large_pitch_roll_error_intervals"]["dt_join"],
    )
    dt = stop - CxoTime(large_pitch_roll_error["datestart"])
    large_pitch_roll_error["recent"] = dt < args.highlight_recent_days * u.day

    pitch_roll_time_outfiles = []
    for interval in large_pitch_roll_error[:CONFIG["max_pitch_roll_time_plots"]]:
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

    # Write out the HTML report
    template = jinja2.Template((data_dir / "index.html").read_text())
    context = {
        "date_last": date_last.date[:-4],
        "date_last_local": date_last_local_fmt,
        "stop": stop.date[:-4],
        "days_long_term": args.days_long_term,
        "days_recent": args.days_recent,
        "large_pitch_roll_error": large_pitch_roll_error,
        "config": CONFIG,
        "axes": ["pitch", "roll"],
        "pitch_roll_time_outfiles": pitch_roll_time_outfiles,
    }
    txt = template.render(**context)
    (outdir / "index.html").write_text(txt)

    # Copy the historical plot files
    outdir_hist = outdir / "fss_prim_hist"
    outdir_hist.mkdir(parents=True, exist_ok=True)
    for inpath in (data_dir / "fss_prim_hist").glob("*"):
        outpath = outdir_hist / inpath.name
        if not outpath.exists():
            shutil.copy2(inpath, outpath)


def get_formatted_local_time(date: CxoTimeLike):
    date = CxoTime(date)
    date_last_local = date.datetime.replace(tzinfo=timezone.utc).astimezone(tz=None)
    date_last_local_fmt = date_last_local.strftime("%Y %a %b %d %I:%M:%S %p %Z")
    return date_last_local_fmt


if __name__ == "__main__":
    main()
