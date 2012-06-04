import asciitable
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

from bad_times import bad_times

SAFEMODE_2012150 = '2012:150:03:33:29'

plt.rc('legend', fontsize=10)

def plot_2008246_event(id='a', savefig=False):
    if id == 'a':
        start = '2008:246:02:00:00'
        stop = '2008:246:03:00:00'
    else:
        start = '2008:246:13:20:00'
        stop = '2008:246:13:45:00'
    dat = fetch.MSIDset(['aoalpang', 'aosunprs', 'pitch'],
                        start, stop)
    dat.interpolate(1.025, filter_bad=True)
    plt.figure(6)
    plt.clf()
    plt.subplot(2, 1, 1)
    plot_cxctime(dat['aoalpang'].times, dat['aoalpang'].vals, '.', label='Sun not present')
    plot_cxctime(dat['aoalpang'].times, dat['aoalpang'].vals)
    ok = dat['aosunprs'].vals == 'SUN '
    plot_cxctime(dat['aoalpang'].times[ok], dat['aoalpang'].vals[ok],
                 '.r', ms=12, mec='r', label='Sun present')
    plt.title('Bad Alpha angles with sun presence on 2008:246')
    plt.ylabel('AOALPANG (deg)')
    plt.legend(loc='upper left')
    plt.grid()
    plt.subplot(2, 1, 2)
    plot_cxctime(dat['pitch'].times, dat['pitch'].vals, '.')
    plt.title('Pitch angle')
    plt.ylabel('pitch (deg)')
    plt.grid()
    plt.tight_layout()
    if savefig:
        plt.savefig('event_2008246{}.png'.format(id))


def plot_altitude(out, ephems=None, angle_err_lim=8.0, savefig=False):
    times = out['times']
    alpha_err = out['alpha'] - out['roll']
    alpha_sun = out['alpha_sun']
    beta_sun = out['beta_sun']
    bad = alpha_sun & beta_sun & (abs(alpha_err) > angle_err_lim)
    if ephems is None:
        ephems = fetch.Msidset(['orbitephem1_*'], times[0], times[-1])
    ephems.interpolate(out['times'][1] - out['times'][0])
    alt_km = np.sqrt(ephems['orbitephem1_x'].vals ** 2 +
                     ephems['orbitephem1_y'].vals ** 2 +
                     ephems['orbitephem1_z'].vals ** 2) / 1000.0
    idxs = np.searchsorted(ephems.times, times[bad])
    idxs = idxs[idxs < len(ephems.times)]
    plt.figure(10)
    plt.clf()
    plot_cxctime(ephems.times, alt_km, ',', color='c', mec='c')
    plot_cxctime(ephems.times[idxs], alt_km[idxs], '.',
                color='b', mec='b', ms=3)
    plt.title('Orbit radius for bad FSS data')
    plt.ylabel('Orbit radius (km)')
    if savefig:
        plt.savefig('orbit_bad_fss.png')
    return ephems, alt_km, bad, idxs


def plot_fss_temp(out, tfss=None, angle_err_lim=8.0, savefig=False):
    times = out['times']
    alpha_err = out['alpha'] - out['roll']
    alpha_sun = out['alpha_sun']
    beta_sun = out['beta_sun']
    bad = alpha_sun & beta_sun & (abs(alpha_err) > angle_err_lim)
    if tfss is None:
        tfss = fetch.Msidset(['tfssbkt1'], times[0], times[-1])
        tfss.interpolate(out['times'][1] - out['times'][0])
    idxs = np.searchsorted(tfss.times, times[bad])
    idxs = idxs[idxs < len(tfss.times)]
    plt.figure(11)
    plt.clf()
    plot_cxctime(tfss.times, tfss['tfssbkt1'].vals, ',', color='c', mec='c')
    plot_cxctime(tfss.times[idxs], tfss['tfssbkt1'].vals[idxs], '.',
                color='b', mec='b', ms=3)
    plt.title('FSS temperature for bad FSS data')
    plt.ylabel('Temperature (degF)')
    if savefig:
        plt.savefig('fss_temp_bad_fss.png')
    return tfss


def plot_pitches(out, angle_err_lim=8.0, savefig=False):
    times= out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    alpha_sun = out['alpha_sun']
    beta_err = out['beta'] - out['pitch']
    beta_sun = out['beta_sun']

    for i, title, xlabel, ylabel in (
        (1, 'Pitch for bad value & sun presence True', None, 'Pitch (deg)'),
        (2, 'Pitch when alpha sun presence is False', None, 'Pitch (deg)'),
        (3, 'Pitch when beta sun presence is False', None, 'Pitch (deg)')):
        plt.figure(i)
        plt.clf()
        plt.grid()
        plt.title(title)
        plt.ylabel(ylabel)
        if xlabel:
            plt.xlabel(xlabel)

    zipvals = zip((~out['kalman'], out['kalman']),
                  ('c', 'r'),
                  ('b', 'r'),
                  ('Not Kalman (cyan)', 'Kalman (red)'))
    for filt, col1, col2, label in zipvals:
        plt.figure(1)
        ok = filt & ~alpha_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1, label=label)

        plt.figure(2)
        ok = filt & ~alpha_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1, label=label)

        plt.figure(3)
        ok = filt & ~beta_sun
        plot_cxctime(times[ok], pitch[ok], ',', color=col1, mec=col1, label=label)

        plt.figure(1)
        ok = filt & alpha_sun & beta_sun & (abs(alpha_err) > angle_err_lim)
        if sum(ok) > 0:
            plot_cxctime(times[ok], pitch[ok], 'o', color=col2, mec=col2,
                         ms=3, label='Bad & sun presence True')

    
    figure(1)
    plot_cxctime([DateTime('2012:150:03:33:00').secs], [139.1], 'x', color='r', mec='r',
                 ms=7, mew=2, label="Safe mode 2012:150")

    suffs = ('bad_alpha_sun', 'alpha_no_sun', 'beta_no_sun')
    for i, suff in enumerate(suffs):
        plt.figure(i + 1)
        x0, x1 = plt.xlim()
        dx = (x1 - x0) / 20
        plt.xlim(x0 - dx, x1 + dx)
        plt.ylim(133.5, 144.5)
        plt.legend(loc='lower right')
        if savefig:
            plt.savefig('pitch_' + suff + '.png')


def plot_angle_err(out, axis='alpha'):
    taxis = axis.title()
    sc = {'alpha': 'roll', 'beta': 'pitch'}
    ok = out['alpha_sun'] & out['beta_sun']
    nok = ~ok
    plt.grid()
    plt.title('')
    plot(out['pitch'][nok], out[axis][nok] - out[sc[axis]][nok], ',b', mec='b')
    if axis == 'beta':
        plt.xlabel('Pitch (deg)')
    plt.ylabel('Angle err (deg)')
    plt.title('{} angle error vs. pitch (AOSUNPRS=NSUN)'.format(taxis))


def plot_angle_errs(out, savefig=False):
    plt.figure(4)
    clf()
    plt.subplot(2, 1, 1)
    plot_angle_err(out, axis='alpha')
    plt.subplot(2, 1, 2)
    plot_angle_err(out, axis='beta')
    if savefig:
        plt.savefig('angle_err.png')


def get_data(start='2005:001', stop=SAFEMODE_2012150, interp=32.8,
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

    out['times'][:] = x['pitch'].times[ok]
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

