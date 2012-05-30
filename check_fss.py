from itertools import count

import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime
import Ska.engarchive.fetch_eng as fetch

from bad_times import bad_times

def make_plots(out, alpha_err_lim=5.0):
    times= out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    alpha_err1 = abs(alpha_err) > alpha_err_lim

    for i, title, xlabel, ylabel in (
        (1, 'Roll error vs. pitch', 'Pitch (deg)', 'Roll err (deg)'),
        (2, 'Roll error vs. time', None, 'Roll err (deg)'),
        (3, 'Glitch pitch vs. time', None, 'Pitch (deg)')):
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
        aerr = alpha_err1 & filt
        sperr = out['aosunprs'] & aerr
        n_sperr = sum(sperr)
        plt.plot(out['pitch'][aerr], alpha_err[aerr], ',', color=col1, mec=col1)
        plt.plot(out['pitch'][::20], alpha_err[::20], ',', color=col1, mec=col1)
        if n_sperr:
            plt.plot(out['pitch'][sperr], alpha_err[sperr], '.',
                     color=col2, mec=col2)

        figure(2)
        plot_cxctime(times[aerr], alpha_err[aerr], ',', color=col1, mec=col1)
        if n_sperr:
            plot_cxctime(times[sperr], alpha_err[sperr], '.',
                         color=col2, mec=col2)

        figure(3)
        plot_cxctime(times[aerr], pitch[aerr], ',', color=col1, mec=col1)
        if n_sperr:
            plot_cxctime(times[sperr], pitch[sperr], '.r',
                         color=col2, mec=col2)

    for i in range(1, 4):
        plt.figure(i)
        x0, x1 = plt.xlim()
        dx = (x1 - x0) / 20
        plt.xlim(x0 - dx, x1 + dx)

def get_data(start='2005:001', stop='2012:144', interp=32.8,
             pitch0=134, pitch1=144):
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aosunprs')
    print 'fetching data'
    x = fetch.MSIDset(msids, start, stop)
    print 'starting interpolate'
    x.interpolate(interp, filter_bad=False)

    for msid in x.values():
        print 'filter_bad_times', msid.msid
        msid.filter_bad_times(table=bad_times)

    ok = ((x['pitch'].vals > pitch0) &
          (x['pitch'].vals < pitch1))

    nvals = np.sum(ok)
    bads = np.zeros(nvals, dtype=bool)
    for msid in x.values():
        print msid.msid, sum(msid.bads[ok])
        bads = bads | msid.bads[ok]
    ok[ok] = ok[ok] & ~bads

    nvals = np.sum(ok)
    colnames = ('times', 'pitch', 'roll', 'alpha', 'beta',
                'aosunprs', 'aopssupm', 'kalman')
    dtypes = ('f8', 'f4', 'f4', 'f4', 'f4', 'bool', 'bool', 'bool')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))

    out['times'][:] = x.times[ok]
    out['pitch'][:] = x['pitch'].vals[ok]
    out['roll'][:] = x['roll'].vals[ok]
    out['alpha'][:] = -x['aoalpang'].vals[ok]
    out['beta'][:] = 90 - x['aobetang'].vals[ok]
    out['aosunprs'][:] = x['aosunprs'].vals[ok] == 'SUN '
    out['aopssupm'][:] = x['aopssupm'].vals[ok] == 'ACT '
    out['kalman'][:] = ((x['aoacaseq'].vals[ok] == 'KALM') &
                        (x['aopcadmd'].vals[ok] == 'NPNT'))
    return out

# np.save('fss4.npy', out)

# pitch_err = pitch - pitch_fss
# roll_err = roll - roll_fss

# figure(1)
# clf()
# plot(pitch, pitch_err, ',')
# figure(2)
# clf()
# plot(roll, roll_err, ',')
# figure(3)
# clf()
# plot(pitch,roll_err, ',')
# plot(pitch[sme], roll_err[sme], 'or', mec='r')

# bad = abs(roll_err) > 5
# figure(4)
# clf()
# plot_cxctime(times[bad], roll_err[bad], ',')

# figure(5)
# clf()
# scatter(pitch, roll_err, c=times, edgecolor='none')
