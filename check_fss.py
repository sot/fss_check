from itertools import count

import asciitable
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

from bad_times import bad_times

def make_plots(out, angle_err_lim=8.0):
    times= out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    alpha_sun = out['alpha_sun']
    beta_err = out['beta'] - out['pitch']
    beta_sun = out['beta_sun']

    for i, title, xlabel, ylabel in (
        (1, 'Pitch vs. time for bad alpha values', None, 'Pitch (deg)'),
        (2, 'Pitch (alpha no sun) vs. time', None, 'Pitch (deg)'),
        (3, 'Pitch (beta no sun) vs. time', None, 'Pitch (deg)')):
        plt.figure(i)
        plt.clf()
        plt.grid()
        plt.title(title)
        plt.ylabel(ylabel)
        if xlabel:
            plt.xlabel(xlabel)

    zipvals = zip((~out['kalman'], out['kalman']),
                  ('c', 'r'),
                  ('b', 'r'))
    for filt, col1, col2 in zipvals:
        plt.figure(1)
        ok = filt & ~alpha_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1)

        continue
        plt.figure(2)
        ok = filt & alpha_sun & (abs(alpha_err) > angle_err_lim)
        if sum(ok) > 0:
            plot_cxctime(times[ok], pitch[ok], ',', color='k', mec='k')
        ok = filt & ~alpha_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1)

        plt.figure(3)
        ok = filt & beta_sun & (abs(beta_err) > angle_err_lim)
        if sum(ok) > 0:
            plot_cxctime(times[ok], pitch[ok], ',', color='k', mec='k')
        ok = filt & ~beta_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1)

    plt.figure(1)
    ok = alpha_sun & (abs(alpha_err) > angle_err_lim)
    if sum(ok) > 0:
        plt.set_cmap(plt.jet())
        plotdates = cxctime2plotdate(times[ok])
        plt.scatter(plotdates, pitch[ok], marker='o',
                    c=np.abs(alpha_err[ok]),
                    edgecolor='none',
                    s=80,
                    vmin=5, vmax=40)
        plt.colorbar()

    for i in range(2, 4):
        plt.figure(i)
        x0, x1 = plt.xlim()
        dx = (x1 - x0) / 20
        plt.xlim(x0 - dx, x1 + dx)
        plt.ylim(133.5, 144.5)

def get_data(start='2005:001', stop='2012:144', interp=32.8,
             pitch0=134, pitch1=144):
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aoalpsun', 'aobetsun')
    print 'fetching data'
    x = fetch.MSIDset(msids, start, stop)

    # Resample MSIDset (values and bad flags) onto a common time sampling
    print 'starting interpolate'
    x.interpolate(interp, filter_bad=False)

    # Remove data during times of known bad or anomalous data
    for msid in x.values():
        print 'filter_bad_times', msid.msid
        filter_bad_times(msid, table=bad_times)

    # Select data only in a limited pitch range
    ok = ((x['pitch'].vals > pitch0) &
          (x['pitch'].vals < pitch1))

    # Determine the logical-or of bad values for all MSIDs and use this
    # to further filter the data sample
    nvals = np.sum(ok)
    bads = np.zeros(nvals, dtype=bool)
    for msid in x.values():
        # Ignore sun position monitor for bad data because it is frequently
        # bad (not available in certain subformats including SSR)
        if msid.MSID == 'AOPSSUPM':
            continue
        print msid.msid, sum(msid.bads[ok])
        bads = bads | msid.bads[ok]
    ok[ok] = ok[ok] & ~bads

    nvals = np.sum(ok)
    colnames = ('times',
                'pitch', 'roll', 'alpha', 'beta',
                'alpha_sun', 'beta_sun', 'spm_act', 'spm_act_bad', 'kalman')
    dtypes = ('f8',
              'f4', 'f4', 'f4', 'f4',
              'bool', 'bool', 'bool', 'bool', 'bool', 'bool')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))

    out['times'][:] = x.times[ok]
    out['pitch'][:] = x['pitch'].vals[ok]
    out['roll'][:] = x['roll'].vals[ok]
    out['alpha'][:] = -x['aoalpang'].vals[ok]
    out['beta'][:] = 90 - x['aobetang'].vals[ok]
    out['alpha_sun'][:] = x['aoalpsun'].vals[ok] == 'SUN '
    out['beta_sun'][:] = x['aobetsun'].vals[ok] == 'SUN '
    out['spm_act'][:] = x['aopssupm'].vals[ok] == 'ACT '
    out['spm_act_bad'][:] = x['aopssupm'].bads[ok]
    out['kalman'][:] = ((x['aoacaseq'].vals[ok] == 'KALM') &
                        (x['aopcadmd'].vals[ok] == 'NPNT'))
    return out

def filter_bad_times(msid_self, start=None, stop=None, table=None):
    """Filter out intervals of bad data in the MSID object.

    There are three usage options:

    - Supply no arguments.  This will use the global list of bad times read
      in with fetch.read_bad_times().
    - Supply both ``start`` and ``stop`` values where each is a single
      value in a valid DateTime format.
    - Supply an ``table`` parameter in the form of a 2-column table of
      start and stop dates (space-delimited) or the name of a file with
      data in the same format.

    The ``table`` parameter must be supplied as a table or the name of a
    table file, for example::

      bad_times = ['2008:292:00:00:00 2008:297:00:00:00',
                   '2008:305:00:12:00 2008:305:00:12:03',
                   '2010:101:00:01:12 2010:101:00:01:25']
      msid.filter_bad_times(table=bad_times)
      msid.filter_bad_times(table='msid_bad_times.dat')

    :param start: Start of time interval to exclude (any DateTime format)
    :param stop: End of time interval to exclude (any DateTime format)
    :param table: Two-column table (start, stop) of bad time intervals
    """
    if table is not None:
        bad_times = asciitable.read(table, Reader=asciitable.NoHeader,
                                    names=['start', 'stop'])
    elif start is None and stop is None:
        bad_times = msid_bad_times.get(msid_self.MSID, [])
    elif start is None or stop is None:
        raise ValueError('filter_times requires either 2 args '
                         '(start, stop) or no args')
    else:
        bad_times = [(start, stop)]

    ok = np.ones(len(msid_self.times), dtype=bool)
    for start, stop in bad_times:
        tstart = DateTime(start).secs
        tstop = DateTime(stop).secs
        if tstart > tstop:
            raise ValueError("Start time %s must be less than stop time %s"
                             % (start, stop))

        if tstop < msid_self.times[0] or tstart > msid_self.times[-1]:
            continue

        i0, i1 = np.searchsorted(msid_self.times, [tstart, tstop])
        ok[i0:i1 + 1] = False

    colnames = (x for x in msid_self.colnames)
    for colname in colnames:
        attr = getattr(msid_self, colname)
        if isinstance(attr, np.ndarray):
            setattr(msid_self, colname, attr[ok])

