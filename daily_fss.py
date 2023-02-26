#!/usr/bin/env python

import argparse
import os
import sys

import matplotlib

matplotlib.use('Agg')  # noqa

from cheta import fetch  # noqa: E402
import Ska.File
from Chandra.Time import DateTime

sys.path.insert(0, os.path.dirname(__file__))  # noqa
from check_fss import get_fss_prim_data, plot_pitches, plot_pitches_any_kalman

parser = argparse.ArgumentParser(description='Daily FSS monitor')
parser.add_argument('--out',
                    type=str,
                    default='.',
                    help='Output data directory')
parser.add_argument('--start',
                    type=str,
                    help='Start date')
parser.add_argument('--stop',
                    type=str,
                    help='Stop date')
parser.add_argument('--interp',
                    type=float,
                    default=4.1,
                    help='Telemetry interpolation (secs)')
parser.add_argument('--use-maude',
                    action='store_true',
                    default=False,
                    help='Use MAUDE to get most recent data')
args = parser.parse_args()

if args.use_maude:
    fetch.data_source.set('cxc', 'maude')

start = DateTime(args.start) if args.start else (DateTime() - 360.0)
start_hist = DateTime(args.start) if args.start else DateTime('2010:001')
stop = DateTime(args.stop)

for fss_dir, get_data, start_t0 in (('fss_prim', get_fss_prim_data, start),
                                    # ('fss_sec', get_fss_sec_data, start),
                                    # ('fss_prim_hist', get_fss_prim_data, start_hist)
                                    ):
    primary = 'prim' in fss_dir
    start = DateTime(args.start or start_t0)
    print('Processing', fss_dir, start.date, stop.date)
    dat = get_data(start, stop, interp=args.interp)
    with Ska.File.chdir(os.path.join(args.out, fss_dir)):
        print(' plot_pitches')
        plot_pitches(dat, savefig=True, start=start, stop=stop, primary=primary,
                     start_suffix=args.start)

    start = DateTime(args.start or start_t0)
    print('Processing', fss_dir, start.date, stop.date)
    dat = get_data(start, stop, interp=args.interp)
    with Ska.File.chdir(os.path.join(args.out, fss_dir)):
        print(' plot_pitches_any_kalman')
        plot_pitches_any_kalman(dat, savefig=True, start=start, stop=stop,
                                primary=primary, start_suffix=args.start)
