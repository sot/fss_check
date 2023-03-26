# Licensed under a 3-clause BSD style license - see LICENSE.rst

import argparse
import os
from pathlib import Path

import matplotlib
import Ska.File
from Chandra.Time import DateTime
from fss_check import get_fss_prim_data, plot_pitches, plot_pitches_any_kalman
from cheta import fetch_eng as fetch


def get_parser():
    parser = argparse.ArgumentParser(
        description="Daily Fine Sun Sensor monitoring and trending"
    )
    parser.add_argument("--out", type=str, default=".", help="Output data directory")
    parser.add_argument("--start", type=str, help="Start date")
    parser.add_argument("--stop", type=str, help="Stop date")
    parser.add_argument(
        "--interp", type=float, default=4.1, help="Telemetry interpolation (secs)"
    )
    parser.add_argument(
        "--use-maude",
        action="store_true",
        default=False,
        help="Use MAUDE to get most recent data",
    )
    return parser


def main(args=None):
    matplotlib.use("Agg")  # noqa

    args = get_parser().parse_args(args)

    if args.use_maude:
        fetch.data_source.set("cxc", "maude allow_subset=False")

    start = DateTime(args.start) if args.start else (DateTime() - 7.0)
    stop = DateTime(args.stop)
    print(f"Processing {start.date} to {stop.date}")

    start = DateTime(args.start)
    print("Processing fss_prim", start.date, stop.date)
    dat = get_fss_prim_data(start, stop, interp=args.interp, use_maude=args.use_maude)

    outdir = Path(args.out) / "fss_prim"
    outdir.mkdir(parents=True, exist_ok=True)

    with Ska.File.chdir(outdir):
        print(" plot_pitches")
        plot_pitches(
            dat,
            savefig=True,
            start=start,
            stop=stop,
            primary=True,
            start_suffix=args.start,
        )

        print(" plot_pitches_any_kalman")
        plot_pitches_any_kalman(
            dat,
            savefig=True,
            start=start,
            stop=stop,
            primary=True,
            start_suffix=args.start,
        )


if __name__ == "__main__":
    main()
