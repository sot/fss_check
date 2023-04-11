# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging
from typing import Optional

import kadi.events.query  # noqa: F401 # For kadi.events.models import (django setup)
import matplotlib.pyplot as plt
import matplotlib.style
import numpy as np
from astropy.table import Table
from cheta.utils import logical_intervals
from cxotime import CxoTime, CxoTimeLike
from kadi import events
from kadi.events.models import fuzz_states
from Ska.Matplotlib import cxctime2plotdate, plot_cxctime, set_min_axis_range

from fss_check.fss_utils import CONFIG, get_spm_pitch_roll

matplotlib.style.use("bmh")

plt.rc("legend", fontsize=10)
plt.rc("lines", markersize=4.0)

events.eclipses.pad_interval = 1000

logger = logging.getLogger("fss_check")


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


def plot_pitch_for_data_with_large_errors(
    dat: Table,
    start: Optional[CxoTimeLike] = None,
    stop: Optional[CxoTimeLike] = None,
    outfile: Optional[str] = None,
    axis: str = "roll",
):
    """Plot pitch for all points where `axis` value error > angle_err_lim.

    Cyan points are with no sun presence, red are with sun presence.
    Unlike plot_pitches() below there is no distinction made based
    on the kalman state.
    """
    logger.info("Running plot_pitch_for_data_with_large_errors")
    times = dat["times"]
    pitch = dat["pitch"]
    pitch_roll_err = dat[f"{axis}_err"]
    sun = dat["alpha_sun"] & dat["beta_sun"]

    vals = [
        (~sun, "c.", "c", 1.0, "No Sun Presense (error > 8.0 deg)", 8.0),
        (sun, "bo", "b", 0.5, "Sun Presense (2.0 < error <= 4.0 deg)", 2.0),
        (sun, "mo", "m", 0.7, "Sun Presense (4.0 < error <= 8.0 deg)", 4.0),
        (sun, "ro", "r", 1.0, "Sun Presense (error > 8.0 deg)", 8.0),
    ]
    plt.figure()
    for filt, mark, mec, alpha, label, err_min in vals:
        ok = (abs(pitch_roll_err) > err_min) & filt
        if np.any(ok):
            logger.info(f" {label} {np.count_nonzero(ok)}")
            plot_cxctime(times[ok], pitch[ok], mark, mec=mec, label=label, alpha=alpha)

    plt.legend(loc="upper left", fancybox=True, framealpha=0.9, fontsize="small")
    plt.title(f"Pitch for {axis} error > threshold")
    plt.ylabel("Pitch (deg)")

    set_plot_limits(start, stop)

    if outfile is not None:
        plt.savefig(outfile)


def get_large_pitch_roll_error_intervals(
    dat: Table,
    axis: str = "roll",
    max_pitch: float = 135,
    max_err: float = 2.0,
    dt_join: float = 100,
) -> Table:
    """Return Table of points with large pitch or roll errors.

    :param dat: data table
    :param axis: 'roll' or 'pitch'
    :param max_pitch: max pitch value (deg)
    :param max_err: max pitch or roll error (deg)
    :param dt_join: time delta (secs) to join intervals
    :returns: table of intervals with large pitch or roll errors
    """
    sun = dat["alpha_sun"] & dat["beta_sun"]
    ok = (np.abs(dat[f"{axis}_err"]) > max_err) & sun & (dat["pitch"] < max_pitch)

    intervals = logical_intervals(dat["times"], ok, max_gap=33)
    intervals = fuzz_states(intervals, dt_join)
    # intervals["duration"].format = ".1f"
    # intervals["tstart"].format = ".1f"
    # intervals["tstop"].format = ".1f"

    return intervals


def plot_pitches(
    dat,
    angle_err_lim=8.0,
    savefig=False,
    start=None,
    stop=None,
    primary=True,
    start_suffix="",
):
    """LEGACY: Plot pitch for all points where roll error > angle_err_lim.

    This is a legacy function that is not used in the current fss_check.
    """
    times = dat["times"]
    pitch = dat["pitch"]
    roll_err = dat["roll_fss"] - dat["roll"]
    alpha_sun = dat["alpha_sun"]
    beta_sun = dat["beta_sun"]

    for i, title, xlabel, ylabel in (
        (
            1,
            "Pitch for alpha error > {} deg".format(angle_err_lim),
            None,
            "Pitch (deg)",
        ),
        (2, "Pitch when alpha sun presence is False", None, "Pitch (deg)"),
        (3, "Pitch when beta sun presence is False", None, "Pitch (deg)"),
    ):
        plt.figure(i)
        plt.clf()
        plt.grid()
        plt.title(title)
        plt.ylabel(ylabel)
        if xlabel:
            plt.xlabel(xlabel)

    zipvals = zip(
        (~dat["kalman"], dat["kalman"]),
        (
            dict(color="c", mec="c"),  # Not Kalman, No sun presence
            dict(color="r", mec="r"),
        ),  # Kalman, No sun presence
        (
            dict(color="b", mec="b", fmt="o"),  # Not Kalman, Sun presence
            dict(color="r", mec="r", fmt="x", mew=2),
        ),  # Kalman, Sun presence
        ("Not Kalman (cyan)", "Kalman (red)"),
    )
    sun_presence = alpha_sun & beta_sun
    bad_value = abs(roll_err) > angle_err_lim
    for filt, opt1, opt2, label in zipvals:
        plt.figure(1)
        ok = filt & bad_value
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ",", label=label, **opt1)

        ok = filt & sun_presence & bad_value
        if np.any(ok):
            plot_cxctime(
                times[ok], pitch[ok], label=label + " & Sun Presence True", **opt2
            )

        plt.figure(2)
        ok = filt & ~alpha_sun
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ",", label=label, **opt1)

        plt.figure(3)
        ok = filt & ~beta_sun
        if np.any(ok):
            plot_cxctime(times[ok], pitch[ok], ",", label=label, **opt1)

    suffs = ("bad_alpha_sun", "alpha_no_sun", "beta_no_sun")
    for i, suff in enumerate(suffs):
        plt.figure(i + 1)

        set_plot_limits(start, stop)

        plt.legend(loc="best")

        if savefig:
            ident = savefig if isinstance(savefig, str) else ""
            plt.savefig("pitch_" + ident + suff + start_suffix + ".png")


def plot_delta_vs_pitch_roll(dat: Table, outfile: Optional[str] = None):
    """Plot delta pitch and delta roll vs. pitch and roll.

    This results in a 2x2 plot with the following subplots:
    - delta pitch vs. pitch
    - delta roll vs. pitch
    - delta pitch vs. roll
    - delta roll vs. roll

    Parameters
    ----------
    dat : astropy.table.Table
        Table of FSS data
    outfile : str, optional
        Filename to write plot to.
    """
    logger.info("Running plot_delta_vs_pitch_roll")
    marker = "."
    marker_size = 1

    ok_pitch = dat["pitch"] < 135
    _, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(8, 6))

    def plot_one(ax, xaxis, yaxis):
        channel = {"pitch": "beta", "roll": "alpha"}[yaxis]
        ok = ok_pitch & dat[channel + "_sun"]

        pitch = dat["pitch"][ok]
        roll = dat["roll"][ok]
        pitch_fss = dat["pitch_fss"][ok]
        roll_fss = dat["roll_fss"][ok]
        d_pitch = (pitch_fss - pitch).clip(-3, 3)
        d_roll = (roll_fss - roll).clip(-3, 3)

        x = pitch if xaxis == "pitch" else roll
        y = d_pitch if yaxis == "pitch" else d_roll

        for mask, color, label in zip(
            (~dat["kalman"][ok], dat["kalman"][ok]),
            ("C0", "C1"),
            ("Kalman", "Not Kalman"),
        ):
            ax.plot(x[mask], y[mask], marker, color=color, ms=marker_size, label=label)
            bad = np.abs(y[mask]) > 1.5
            if np.any(bad):
                ax.plot(x[mask][bad], y[mask][bad], marker, color=color)
        ax.margins(0.05)
        ax.set_xlabel(f"{xaxis.capitalize()} (deg)")
        ax.set_ylabel(f"Delta {yaxis} (deg)")
        ax.set_title(f"FSS {yaxis} - OBC {yaxis} vs. {xaxis}")
        x0, x1 = ax.get_xlim()
        ax.hlines(
            [1.5, -1.5], xmin=x0, xmax=x1, color="r", linestyle="dashed", alpha=0.5
        )
        ax.set_xlim(x0, x1)
        if xaxis == "pitch" and yaxis == "roll":
            ax.legend(loc="upper left")

    # fmt: off
    plot_one(ax0, "roll", "pitch")
    plot_one(ax1, "pitch", "pitch")
    plot_one(ax2, "roll", "roll")
    plot_one(ax3, "pitch", "roll")
    # fmt: on

    year = CxoTime(dat["times"][0]).date[:4]
    plt.suptitle(f"FSS data for {year} through {CxoTime(dat['times'][-1]).date}")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile)


def plot_roll_pitch_vs_time(dat, start, stop):
    """Plot roll and pitch (OBC and FSS) vs. time.

    This is used to make detail plots of short time intervals.

    Parameters
    ----------
    dat : astropy.table.Table
        Table of FSS data
    start : CxoTimeLike
        Start time
    stop : CxoTimeLike
        Stop time
    """
    start = CxoTime(start)
    stop = CxoTime(stop)
    i0, i1 = np.searchsorted(dat["times"], [start.secs, stop.secs])
    dat = dat[i0:i1]

    sun_present = dat["alpha_sun"] & dat["beta_sun"] & (dat["pitch"] < 135)
    roll_fss = dat["roll_fss"].copy()
    roll_err = dat["roll_fss"] - dat["roll"]
    roll_fss[~sun_present] = np.nan
    roll_err[~sun_present] = np.nan
    pitch_fss = dat["pitch_fss"].copy()
    pitch_err = dat["pitch_fss"] - dat["pitch"]
    pitch_fss[~sun_present] = np.nan
    pitch_err[~sun_present] = np.nan

    _, axs = plt.subplots(2, 2, figsize=(10, 4), sharex=True)
    plot_cxctime(
        dat["times"],
        roll_fss,
        ".-",
        ax=axs[0, 0],
        lw=0.25,
        color="C1",
        label="Roll FSS",
    )
    plot_cxctime(
        dat["times"], dat["roll"], "-", ax=axs[0, 0], color="C0", label="Roll OBC"
    )
    plot_cxctime(
        dat["times"],
        roll_err,
        ".-",
        color="C1",
        ax=axs[1, 0],
        ms=3,
        lw=0.25,
        label="Roll error",
    )
    plot_cxctime(
        dat["times"],
        pitch_fss,
        ".-",
        ax=axs[0, 1],
        lw=0.5,
        color="C1",
        label="Pitch FSS",
    )
    plot_cxctime(
        dat["times"], dat["pitch"], "-", ax=axs[0, 1], color="C0", label="Pitch OBC"
    )
    plot_cxctime(
        dat["times"],
        pitch_err,
        ".",
        color="C1",
        ax=axs[1, 1],
        ms=3,
        label="Pitch error",
    )

    axs[0, 0].legend(loc="best")
    axs[0, 1].legend(loc="best")
    axs[0, 0].set_ylabel("Roll (deg)")
    axs[1, 0].set_ylabel("Roll error (deg)")
    axs[0, 1].set_ylabel("Pitch (deg)")
    axs[1, 1].set_ylabel("Pitch error (deg)")
    set_min_axis_range(axs[1, 0], 1.0)
    set_min_axis_range(axs[1, 1], 1.0)
    plt.tight_layout()


def plot_pitch_roll_spm_mp_constraints(dat):
    from ska_sun import ROLL_TABLE

    plt.figure(figsize=(12, 8))
    pitch, roll = get_spm_pitch_roll()
    plt.plot(pitch, roll, color="C1", label="SPM limit")
    ok = dat["alpha_sun"] & dat["beta_sun"] & (dat["pitch"] < 135)
    dok = dat[ok]
    plt.plot(dok["pitch_fss"], dok["roll_fss"], ".", ms=1, alpha=0.5, color="C0")
    bad = np.abs(dok["roll_fss"] - dok["roll"]) > 1.0
    plt.plot(
        dok["pitch_fss"][bad],
        dok["roll_fss"][bad],
        ".",
        color="C1",
        ms=4,
        label="FSS roll err > 1 deg",
    )
    bad = np.abs(dok["roll_fss"] - dok["roll"]) > 1.5
    plt.plot(
        dok["pitch_fss"][bad],
        dok["roll_fss"][bad],
        ".",
        color="r",
        ms=8,
        label="FSS roll err > 1.5 deg",
    )
    plt.plot(
        ROLL_TABLE.val["pitch"],
        ROLL_TABLE.val["rolldev"],
        color="C1",
        lw=1,
        label="Planning limit",
    )
    plt.plot(ROLL_TABLE.val["pitch"], -ROLL_TABLE.val["rolldev"], color="C1", lw=1)
    plt.xlim(None, 140)
    plt.xlabel("Pitch (deg)")
    plt.ylabel("Roll (deg)")
    plt.title(
        "Pitch and roll from FSS data 2022 and 2023 through "
        f'{CxoTime(dat["times"][-1]).date}'
    )
    plt.legend(loc="upper left")


def plot_no_sun_presence(dat):
    """Plot pitch for all points where no sun presence is indicated."""
    times = dat["times"]
    pitch = dat["pitch"]
    sun = dat["alpha_sun"] & dat["beta_sun"]

    p135 = pitch < 135
    ok = ~sun & p135
    plot_cxctime(
        times[ok],
        pitch[ok],
        ".",
        color="C1",
        label="No sun presence, pitch < 135",
        alpha=0.5,
    )

    ok = ~sun
    plot_cxctime(
        times[ok],
        pitch[ok],
        ",",
        color="C0",
        label="No sun presence, pitch < 135",
    )
    last_date = CxoTime(times[-1]).date[:-4]

    plt.legend(loc="best", fancybox=True, framealpha=0.5)
    plt.grid("on")
    plt.title(f"No sun presence (not AOALPSUN or not AOBETSUN) (through {last_date})")
    plt.ylabel("Pitch (deg)")
