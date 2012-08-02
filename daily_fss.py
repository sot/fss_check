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
import check_fss

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

stop = args.stop or DateTime()
start = args.start or stop - 180
dat = check_fss.get_data(start, stop, interp=args.interp)
with Ska.File.chdir(args.out):
    check_fss.plot_pitches(dat, savefig=True)

start = args.start or stop - 365
dat = check_fss.get_data(start, stop, interp=args.interp)
with Ska.File.chdir(args.out):
    check_fss.plot_pitches_any_kalman(dat, savefig=True)
