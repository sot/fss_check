import sys
import os

import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime
from kadi import events

sys.path.insert(0, os.path.dirname(__file__))
from bad_times import bad_times

plt.rc('legend', fontsize=10)
events.eclipses.pad_interval = 1000

def plot_swap_line(primary):
    swap_date = DateTime('2013:130:20:00:00')
    swap_x = cxctime2plotdate([swap_date.secs])
    x0, x1 = plt.xlim()
    y0, y1 = plt.ylim()
    plt.plot([swap_x, swap_x], [y0, y1], '--g', lw=2)
    text_y = y1 - (y1 - y0) * 0.08
    dx = (x1 - x0) * 0.05
    label1, label2 = ('FSS-A', 'FSS-B') if primary else ('FSS-B', 'FSS-A')
    plt.text(swap_x - dx, text_y, label1, ha='right',
             bbox=dict(facecolor='yellow'))
    plt.text(swap_x + dx, text_y, label2, ha='left',
             bbox=dict(facecolor='yellow'))


def set_plot_limits(start, stop):
    if start and stop:
        x0, x1 = cxctime2plotdate([start.secs, stop.secs])
        plt.xlim(x0, x1)

    x0, x1 = plt.xlim()
    dx = (x1 - x0) / 20
    plt.xlim(x0 - dx, x1 + dx)
    y0, y1 = plt.ylim()
    y0 = min(y0, 133.5)
    dy = (y1 - y0) / 20
    plt.ylim(y0 - dy, y1 + dy)


def plot_pitches_any_kalman(out, angle_err_lim=8.0, savefig=False, start=None, stop=None,
                            primary=True):
    """Plot pitch for all points where alpha_err > angle_err_lim.
    Cyan points are with no sun presence, red are with sun presence.
    Unlike plot_pitches() below there is no distinction made based
    on the kalman state.
    """
    times = out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    sun = out['alpha_sun'] & out['beta_sun']
    bad = abs(alpha_err) > angle_err_lim

    zipvals = zip((~sun, sun),
                  ('c.', 'r.'),
                  ('c', 'r'),
                  ('No Sun Presence', 'Sun Presence'))
    plt.figure()
    for filt, mark, mec, label in zipvals:
        ok = bad & filt
        if np.any(ok):
            plot_cxctime(times[bad & filt], pitch[ok], mark,
                         mec=mec, label=label)
    plt.legend(loc='lower left')
    plt.grid('on')
    plt.title("Pitch for alpha error > {} deg".format(angle_err_lim))
    plt.ylabel('Pitch (deg)')

    set_plot_limits(start, stop)

    plot_swap_line(primary)

    if savefig:
        plt.savefig('pitch_bad_alpha.png')


def plot_pitches(out, angle_err_lim=8.0, savefig=False, start=None, stop=None,
                 primary=True):
    times = out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    alpha_sun = out['alpha_sun']
    beta_sun = out['beta_sun']

    for i, title, xlabel, ylabel in (
        (1, "Pitch for alpha error > {} deg".format(angle_err_lim), None,
         'Pitch (deg)'),
        (2, 'Pitch when alpha sun presence is False', None, 'Pitch (deg)'),
        (3, 'Pitch when beta sun presence is False', None, 'Pitch (deg)')):
        plt.figure(i)
        plt.clf()
        plt.grid()
        plt.title(title)
        plt.ylabel(ylabel)
        if xlabel:
            plt.xlabel(xlabel)

    zipvals = zip((~out['kalman'],
                    out['kalman']),
                  (dict(color='c', mec='c'),  # Not Kalman, No sun presence
                   dict(color='r', mec='r')),  # Kalman, No sun presence
                  (dict(color='b', mec='b', fmt='o'), # Not Kalman, Sun presence
                   dict(color='r', mec='r', fmt='x', mew=2)), # Kalman, Sun presence
                  ('Not Kalman (cyan)',
                   'Kalman (red)'))
    sun_presence = alpha_sun & beta_sun
    bad_value = abs(alpha_err) > angle_err_lim
    for filt, opt1, opt2, label in zipvals:
        plt.figure(1)
        ok = filt & bad_value
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ',',
                         label=label, **opt1)

        ok = filt & sun_presence & bad_value
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok],
                         label=label + ' & Sun Presence True',
                         **opt2)

        plt.figure(2)
        ok = filt & ~alpha_sun
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ',',
                         label=label, **opt1)

        plt.figure(3)
        ok = filt & ~beta_sun
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ',',
                         label=label, **opt1)

    suffs = ('bad_alpha_sun', 'alpha_no_sun', 'beta_no_sun')
    for i, suff in enumerate(suffs):
        plt.figure(i + 1)

        set_plot_limits(start, stop)

        plt.legend(loc='best')
        plot_swap_line(primary)

        if savefig:
            ident = savefig if isinstance(savefig, basestring) else ''
            plt.savefig('pitch_' + ident + suff + '.png')


def get_fss_prim_data(start='2011:001', stop=DateTime().date, interp=4.1,
                      pitch0=100, pitch1=144):
    """
    Get data for the primary FSS (FSS-A before ~2013:130:20:00:00, FSS-B after)
    """
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aoalpsun', 'aobetsun')
    print 'fetching data'
    x = fetch.MSIDset(msids, start, stop)

    # Resample MSIDset (values and bad flags) onto a common time sampling
    print 'starting interpolate'
    x.interpolate(interp, filter_bad=False)

    # Remove data during times of known bad or anomalous data (works as of
    # Ska.engarchive 0.19.1)
    x.filter_bad_times(table=bad_times)
    for msid in msids:
        x[msid].remove_intervals(events.eclipses | events.safe_suns)

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
        print msid.msid, np.sum(msid.bads[ok])
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


def get_fss_sec_data(start='2012:230', stop=DateTime().date, interp=4.1,
                     pitch0=100, pitch1=144):
    """
    Get data for the secondary FSS (FSS-B before ~2013:130:20:00:00, FSS-A after)
    """
    msids = ('aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aspefsw2a', 'aspefsw4a', 'aspefsw2b', 'aspefsw4b',
             'ccsdsvcd', 'cotlrdsf')
    print 'fetching data'
    if DateTime(start).date < '2012:230':
        start = '2012:230'
    x = fetch.MSIDset(msids, start, stop)

    # Resample MSIDset (values and bad flags) onto a common time sampling
    # defined by the times when the FSS-secondary values are telemetered.
    print 'starting interpolate'
    interpolate_times(x, times=x['aspefsw2b'].times, filter_bad=False)

    # Remove data during times of known bad or anomalous data (works as of
    # Ska.engarchive 0.19.1)
    x.filter_bad_times(table=bad_times)
    for msid in msids:
        x[msid].remove_intervals(events.eclipses | events.safe_suns)

    # Select data only in a limited pitch range
    pitch_range = ((x['pitch'].vals > pitch0) &
                   (x['pitch'].vals < pitch1))

    # Good data are telemetered at minor frames [0, 32, 64, 96]
    # Bogus data at [16, 48, 72, 104]
    good_mf = np.mod(x['ccsdsvcd'].vals, 32) == 0

    # Select data in PCAD diagnostic subformat (otherwise no FSS-B telemetry)
    pcad_sfmt = x['cotlrdsf'].vals == 'PCAD'

    ok = pitch_range & good_mf & pcad_sfmt

    # Determine the logical-or of bad values for all MSIDs and use this
    # to further filter the data sample
    nvals = np.sum(ok)
    bads = np.zeros(nvals, dtype=bool)
    for msid in x.values():
        # Ignore sun position monitor for bad data because it is frequently
        # bad (not available in certain subformats including SSR)
        if msid.MSID == 'AOPSSUPM':
            continue
        print msid.msid, np.sum(msid.bads[ok])
        bads = bads | msid.bads[ok]
    ok[ok] = ok[ok] & ~bads

    nvals = np.sum(ok)
    colnames = ('times',
                'pitch', 'roll', 'alpha', 'beta',
                'alpha_sun', 'beta_sun', 'kalman')
    dtypes = ('f8',
              'f4', 'f4', 'f4', 'f4',
              'bool', 'bool', 'bool', 'bool', 'bool', 'bool')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))

    out['times'][:] = x['pitch'].times[ok]
    out['pitch'][:] = x['pitch'].vals[ok]
    out['roll'][:] = x['roll'].vals[ok]
    out['alpha'][:] = -x['aspefsw2b'].vals[ok]
    out['beta'][:] = 90 - x['aspefsw4b'].vals[ok]
    out['alpha_sun'][:] = x['aspefsw2a'].vals[ok] == 'SUN '
    out['beta_sun'][:] = x['aspefsw4a'].vals[ok] == 'SUN '
    out['kalman'][:] = ((x['aoacaseq'].vals[ok] == 'KALM') &
                        (x['aopcadmd'].vals[ok] == 'NPNT'))
    return out


def interpolate_times(msidset, times, filter_bad=True):
    """Perform nearest-neighbor interpolation of all MSID values in the set
    to a common time sequence.  The values are updated in-place.

    For each MSID in the set the ``times`` attribute is set to the common
    time sequence.  In addition a new attribute ``times0`` is defined that
    stores the nearest neighbor interpolated time, providing the *original*
    timestamps of each new interpolated value for that MSID.

    By default ``filter_bad`` is True and each MSID has bad data filtered
    *before* interpolation so that the nearest neighbor interpolation only
    finds good data.  In this case (or for ``stat`` values which do not
    have bad values), the ``times0`` attribute can be used to diagnose
    whether the interpolation is meaningful.  For instance large numbers of
    identical time stamps in ``times0`` may be a problem.

    If ``filter_bad`` is set to false then data are not filtered.  In this
    case the MSID ``bads`` values are interpolated as well and provide an
    indication if the interpolated values are bad.

    :param msidset: input MSIDset
    :param times: times at which to interpolate the MSIDset (sec)
    :param filter_bad: filter bad values before interpolating
    """
    import Ska.Numpy

    msids = msidset.values()  # MSID objects in the MSIDset

    # Ensure that tstart / tstop is entirely within the range of available
    # data fetched from the archive.
    max_fetch_tstart = max(msid.times[0] for msid in msids)
    min_fetch_tstop = min(msid.times[-1] for msid in msids)
    istart, istop = np.searchsorted(times, [max_fetch_tstart, min_fetch_tstop])

    msidset.times = times[istart:istop]

    for msid in msids:
        if filter_bad:
            msid.filter_bad()
        indexes = Ska.Numpy.interpolate(np.arange(len(msid.times)),
                                        msid.times, msidset.times,
                                        method='nearest', sorted=True)
        for colname in msid.colnames:
            colvals = getattr(msid, colname)
            if colvals is not None:
                setattr(msid, colname, colvals[indexes])

        # Make a new attribute times0 that stores the nearest neighbor
        # interpolated times.  Then set the MSID times to be the common
        # interpolation times.
        msid.times0 = msid.times
        msid.times = msidset.times
