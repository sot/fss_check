#!/usr/bin/env python

import sys
import os
import argparse

import matplotlib
matplotlib.use('Agg')

import Ska.File
from Chandra.Time import DateTime

# Allow this script to fetch if being run in a development Ska
os.environ['ENG_ARCHIVE'] = '/proj/sot/ska/data/eng_archive'

sys.path.insert(0, os.path.dirname(__file__))
from check_fss import (get_fssa_data, get_fssb_data,
                       plot_pitches, plot_pitches_any_kalman)

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
args = parser.parse_args()

stop = DateTime(args.stop)

for fss_dir, get_data, start_t0 in (('fssa', get_fssa_data, '2012:230'),
                                    ('fssb', get_fssb_data, '2012:230'),
                                    ('fssa_hist', get_fssa_data, '2010:001')):
    start = DateTime(args.start or start_t0)
    print 'Processing', fss_dir, start.date, stop.date
    dat = get_data(start, stop, interp=args.interp)
    with Ska.File.chdir(os.path.join(args.out, fss_dir)):
        print ' plot_pitches'
        plot_pitches(dat, savefig=True)

    start = DateTime(args.start or start_t0)
    print 'Processing', fss_dir, start.date, stop.date
    dat = get_data(start, stop, interp=args.interp)
    with Ska.File.chdir(os.path.join(args.out, fss_dir)):
        print ' plot_pitches_any_kalman'
        plot_pitches_any_kalman(dat, savefig=True)
