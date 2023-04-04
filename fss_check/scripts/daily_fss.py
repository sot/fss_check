# Licensed under a 3-clause BSD style license - see LICENSE.rst

import argparse
import shutil
from pathlib import Path

import astropy.units as u
import jinja2
import matplotlib
import numpy as np
import ska_file
from cheta import fetch_eng as fetch
from cxotime import CxoTime
from datetime import timezone

import fss_check
from fss_check.check_fss import plot_pitches_any_kalman
from fss_check.fss_utils import add_pitch_roll_columns, get_fss_prim_data

DATA_DIR = Path(fss_check.__file__).parent / "data"


def get_parser():
    parser = argparse.ArgumentParser(
        description="Daily Fine Sun Sensor monitoring and trending"
    )
    parser.add_argument("--out", type=str, default=".", help="Output data directory")
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
    return parser


def main(args=None):
    matplotlib.use("Agg")

    args = get_parser().parse_args(args)

    if args.days_recent > args.days_long_term:
        raise ValueError("recent days must be <= long term days")

    if args.use_maude:
        fetch.data_source.set("cxc", "maude allow_subset=False")

    stop = CxoTime(args.stop)
    start_long_term = stop - args.days_long_term * u.day
    start_recent = stop - args.days_recent * u.day

    print(f"Processing {start_long_term.date} to {stop.date}")
    dat = get_fss_prim_data(start_long_term, stop)
    dat = add_pitch_roll_columns(dat)

    date_last = CxoTime(dat["times"][-1])
    date_last_local = date_last.datetime.replace(tzinfo=timezone.utc).astimezone(
        tz=None
    )
    date_last_local_fmt = date_last_local.strftime("%Y %a %b %d %I:%M:%S %p %Z")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    for start, suffix in [(start_long_term, "long_term"), (start_recent, "recent")]:
        i0 = np.searchsorted(dat["times"], start.cxcsec)
        dok = dat[i0:]
        for axis in ["pitch", "roll"]:
            with ska_file.chdir(outdir):
                print(" plot_pitches_any_kalman")
                plot_pitches_any_kalman(
                    dat=dok,
                    savefig=True,
                    start=start,
                    stop=stop,
                    suffix=suffix,
                    axis=axis,
                )

    # Write out the HTML report
    template = jinja2.Template((DATA_DIR / "index.html").read_text())
    context = {
        "date_last": date_last.date[:-4],
        "date_last_local": date_last_local_fmt,
        "stop": stop.date[:-4],
        "days_long_term": args.days_long_term,
        "days_recent": args.days_recent,
    }
    txt = template.render(**context)
    (outdir / "index.html").write_text(txt)

    # Copy the historical plot files
    outdir_hist = outdir / "fss_prim_hist"
    outdir_hist.mkdir(parents=True, exist_ok=True)
    for inpath in (DATA_DIR / "fss_prim_hist").glob("*"):
        outpath = outdir_hist / inpath.name
        if not outpath.exists():
            shutil.copy2(inpath, outpath)


if __name__ == "__main__":
    main()
