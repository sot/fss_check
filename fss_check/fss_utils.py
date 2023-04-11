# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""FSS utilities"""
import os
from pathlib import Path

import cheta.fetch_eng as fetch
import numpy as np
import ska_numpy
import yaml
from astropy.table import Table
from cxotime import CxoTime
from kadi import events
from ska_helpers.utils import LazyDict

SKA = Path(os.environ["SKA"])


def _load_config_file():
    """Load fss_check_config.yaml from current dir or $SKA/data/fss_check"""
    for config_dir, from_ska in [(".", False), (SKA / "data" / "fss_check", True)]:
        path = Path(config_dir) / "fss_check_config.yml"
        if path.exists():
            break

    with open(path) as fh:
        config = yaml.safe_load(fh)

    config['config_path'] = str(path.absolute())
    config['config_from_ska'] = from_ska

    return config


CONFIG = LazyDict(load_func=_load_config_file)

# ### FSS counts to angle calibration on OBC
#
# This uses the on-board formula `angle = arctan(c0 + c1 * counts) + c2)`.
#
# The ground TDB calibration is a 9th-order polynomial that approximates this. For
# all this analysis we use the OBC angle since that is what matters for S/C operations.

FSS_OBC_CAL = {
    "a": {
        "alpha": [-1.190313, 0.00007261314, -0.00001571],
        "beta": [-1.190563, 0.00007264188, -0.00005934],
    },
    "b": {
        "alpha": [-1.190163, 0.00007266802, 0.0001902],
        "beta": [-1.19044, 0.00007265815, -0.0001902],
    },
}

# ### Align FSS to ACA (Chandra body) frame

# FSS-A from OFLS characteristics
ODB_FSS_MISALIGN = np.zeros([4, 4], dtype=np.float64)
ODB_FSS_MISALIGN[1, 1] = 9.999990450374580e-01
ODB_FSS_MISALIGN[2, 1] = -5.327615067743422e-07
ODB_FSS_MISALIGN[3, 1] = 1.381999959551952e-03
ODB_FSS_MISALIGN[1, 2] = 0.0
ODB_FSS_MISALIGN[2, 2] = 9.999999256947376e-01
ODB_FSS_MISALIGN[3, 2] = 3.855003493343671e-04
ODB_FSS_MISALIGN[1, 3] = -1.382000062241829e-03
ODB_FSS_MISALIGN[2, 3] = -3.854999811959735e-04
ODB_FSS_MISALIGN[3, 3] = 9.999989707322665e-01

# FSS-B reference PR-391 from W.S. Davis analysis which determined an alignment
# matrix which minimizes FSS pitch and roll errors.
FSS_TO_ACA = {
    "a": ODB_FSS_MISALIGN[1:, 1:],
    "b": [
        [0.99999147, -0.00348329, -0.00221858],
        [0.00347606, 0.99998865, -0.00325693],
        [0.00222990, 0.00324919, 0.99999224],
    ],
}

# ### TDB polynominal calibration coefficients
#
# These are reversed for use in `np.polyval` and they get used to convert from
# FSS angle in the CXC archive to raw counts.

ALPHA_ANG_PC = [
    -50.01873285,
    0.001739300837,
    4.202023792e-08,
    9.617900394e-12,
    -1.570367949e-15,
    1.756637723e-19,
    -1.118788495e-23,
    3.775679639e-28,
    -6.421573874e-33,
    4.369576818e-38,
][::-1]
BETA_ANG_PC = [
    -50.16789839,
    0.001739575081,
    4.202728171e-08,
    9.631457701e-12,
    -1.573063077e-15,
    1.760037996e-19,
    -1.121133178e-23,
    3.784196473e-28,
    -6.437103687e-33,
    4.380850442e-38,
][::-1]

# Generate a curve corresponding to the ground (TDB) calibration which
# gets numerically inverted to compute on-board counts from ground angle.
COUNTS = np.arange(40000)
ALPHA_ANG = np.polyval(ALPHA_ANG_PC, COUNTS)
BETA_ANG = np.polyval(BETA_ANG_PC, COUNTS)

# Variation of what will be in ska_helpers.utils


def cache_file(fn):
    """Cache the np.ndarray result of a function to a file.

    The cache file is in a local directory called "cache" and is named with the
    function name combined with the arguments and keyword arguments.

    The cache file is a pickle file with the function results.

    Parameters
    ----------
    fn : function

    Returns
    -------
    wrapped : function
    """
    import os
    from pathlib import Path

    def wrapped(*args, **kwargs):
        # define a wrapper that will finally call "fn" with all arguments
        # if cache exists -> load it and return its content
        if args:
            args_str = "_" + "_".join([str(arg) for arg in args])
        else:
            args_str = ""
        if kwargs:
            kwargs_str = "_" + "_".join(
                [f"{key}_{value}" for key, value in kwargs.items()]
            )
        else:
            kwargs_str = ""
        cachedir = Path("cache")
        cachedir.mkdir(exist_ok=True)
        cachefile = cachedir / f"{fn.__name__}{args_str}{kwargs_str}.npz"

        if os.path.exists(cachefile):
            with np.load(cachefile) as fh:
                res = fh["arr_0"]
        else:
            # execute the function with all arguments passed
            res = fn(*args, **kwargs)
            # write to cache file
            np.savez_compressed(cachefile, res)

        return res

    return wrapped


def fss_angle_ground_to_obc(angle, channel, fss="b"):
    """
    Convert ``angle`` from ground telemetry (AOALPANG or AOBETANG) to output ``angle``
    (degrees) in the OBC for the given ``channel`` ('alpha' or 'beta') and ``fss``
    ('a' or 'b').

    Parameters
    ----------
    angle : float, ndarray
        TDB (ground) calibration angle in degrees
    channel : {'alpha', 'beta'}
        FSS channel
    fss : {'a', 'b'}
        FSS identifier

    Returns
    -------
    angle : float, ndarray
        OBC angle in degrees
    """
    y = ALPHA_ANG if channel == "alpha" else BETA_ANG
    f_counts = ska_numpy.interpolate(COUNTS, y, angle)
    counts = ska_numpy.interpolate(COUNTS, y, angle, method="nearest")
    if np.any(np.abs(f_counts - counts) > 0.01):
        raise ValueError("Bad calibration inversion")

    c0, c1, c2 = FSS_OBC_CAL[fss][channel]
    angle = np.degrees(np.arctan(c0 + c1 * counts) + c2)

    return angle


def arccos_clip(x):
    """Return arccos(x) where x may be slightly outside the range [-1, 1]."""
    return np.arccos(x.clip(-1, 1))


def fss_coupled_angles(alpha, beta, fss="b"):
    """
    Convert FSS alpha,beta angles to corrected angles using cross-coupling coeffs.

    See 2015 FSS Calibration Report by W.S. Davis.

    Parameters
    ----------
    alpha : float, ndarray
        FSS alpha angle (deg, rad)
    beta : float, ndarray
        FSS beta angle (deg, rad)
    fss : {'a', 'b'}
        FSS identifier

    Returns
    -------
    alpcor : float, ndarray
        Corrected alpha angle (deg, rad)
    betcor : float, ndarray
        Corrected beta angle (deg, rad)
    """
    if fss == "a":
        alpcor = alpha
        betcor = beta
    elif fss == "b":
        # cross coupling correction parameters - FSS-B
        gcb_caa = 1 + 0.00061451
        gcb_cab = 0.00016852
        gcb_cba = 0.00527152
        gcb_cbb = 1 - 0.00014643
        alpcor = (gcb_caa * alpha) + (gcb_cab * beta)
        betcor = (gcb_cba * alpha) + (gcb_cbb * beta)
    else:
        raise ValueError(f"Invalid FSS: {fss}")

    return alpcor, betcor


def get_pitch_roll_fss(alpha, beta, fss="b", obc=True):
    """Convert FSS alpha, beta angle telemetry to Sun pitch and off-nominal roll.

    Parameters
    ----------
    alpha : float or array-like
        FSS alpha angle (AOALPANG from telemetry, degrees)
    beta : float or array-like
        FSS beta angle (AOBETANG from telemetry, degrees)
    fss : {'a', 'b'}
        FSS identifier
    obc : bool
        If True (default) then use OBC calibration.  If False then use ground.

    Returns
    -------
    pitch : float or ndarray
        Sun pitch angle from FSS data (degrees)
    roll : float or ndarray
        Sun off-nominal roll angle from FSS data (degrees)
    """
    alpha = np.asarray(alpha)
    beta = np.asarray(beta)

    if obc:
        # Convert from ground TDB calibration (9th order polynomial) to OBC calibration
        # (np.arctan(c0 + c1 * counts) + c2)).
        alpha = fss_angle_ground_to_obc(alpha, "alpha", fss=fss)
        beta = fss_angle_ground_to_obc(beta, "beta", fss=fss)

    # Apply cross-coupling correction between channels.
    alpha, beta = fss_coupled_angles(alpha, beta, fss=fss)

    # Adapted from FSS calculations written by A. Arvai
    # https://github.com/sot/cheta/blob/847465a4953/cheta/derived/pcad.py#L337

    # Convert to sun vector in FSS frame.
    alpha_rad = np.deg2rad(np.asarray(alpha, dtype=float))
    beta_rad = np.deg2rad(np.asarray(beta, dtype=float))
    sun_fss = np.array([np.tan(beta_rad), np.tan(alpha_rad), -np.ones_like(alpha_rad)])

    # Convert FSS sun vector to ACA frame.
    sun_aca = FSS_TO_ACA[fss] @ sun_fss
    magnitude = np.linalg.norm(sun_aca, axis=0)
    sun_vec_norm = sun_aca / magnitude
    pitch_fss = np.rad2deg(arccos_clip(sun_vec_norm[0]))
    roll_fss = np.rad2deg(np.arctan2(-sun_vec_norm[1], -sun_vec_norm[2]))

    return pitch_fss, roll_fss


def get_fss_prim_data(start, stop=None, offset=0.9, pitch0=40, pitch1=144):
    """
    Get data for the primary FSS (FSS-A before ~2013:130:20:00:00, FSS-B after)

    Returns structured ndarray with columns::

        "times"
        "pitch"
        "roll"
        "alpang"
        "betang"
        "alpha_sun"
        "beta_sun"
        "kalman"

    Parameters
    ----------
    start : CxoTimeLike
        Start time
    stop : CxoTimeLike
        Stop time (default = start of next year)
    offset : float
        Time offset to apply to FSS data (sec)
    pitch0 : float
        Minimum pitch angle for data (deg)
    pitch1 : float
        Maximum pitch angle for data (deg)

    Returns
    -------
    dat : structured ndarray
    """
    if stop is None:
        start_year = int(CxoTime(start).date[:4])
        stop = f"{start_year + 1}:001"

    msids = (
        "aopcadmd",
        "aoacaseq",
        "pitch",
        "roll",  #
        "aoalpang",
        "aobetang",
        "aoalpsun",
        "aobetsun",
    )
    dat = fetch.MSIDset(msids, start, stop)

    # Apply time offsets to FSS data. This empirically reduces the offset between
    # OBC estimated pitch/roll and FSS pitch/roll during maneuvers.
    for name in ("aoalpang", "aobetang", "aoalpsun", "aobetsun"):
        dat[name].times += offset

    # Resample MSIDset (values and bad flags) onto a common time sampling
    dat.interpolate(times=dat["aoalpang"].times[1:-1], filter_bad=False)

    # Remove data during times of known bad or anomalous data.
    dat.filter_bad_times(table=CONFIG["exclude_intervals"])
    events.eclipses.interval_pad = (1000, 1000)
    events.safe_suns.interval_pad = (0, 100000)
    events.normal_suns.interval_pad = (1000, 10000)
    for msid in msids:
        dat[msid].remove_intervals(
            events.eclipses | events.normal_suns  #  | events.safe_suns
        )

    # Select data only in a limited pitch range
    ok = (
        (dat["pitch"].vals > pitch0)
        & (dat["pitch"].vals < pitch1)
        & (np.isin(dat["aopcadmd"].vals, ["NMAN", "NPNT"]))
    )

    # Determine the logical-or of bad values for all MSIDs and use this
    # to further filter the data sample
    nvals = np.sum(ok)
    bads = np.zeros(nvals, dtype=bool)
    for msid in list(dat.values()):
        bads = bads | msid.bads[ok]
    ok[ok] = ok[ok] & ~bads

    nvals = np.sum(ok)
    colnames = (
        "times",
        "pitch",
        "roll",
        "alpang",
        "betang",
        "alpha_sun",
        "beta_sun",
        "kalman",
    )
    dtypes = (
        "f8",
        "f4",
        "f4",
        "f4",
        "f4",
        "bool",
        "bool",
        "bool",
    )
    # NOTE: leave this output as a structured array so that the file caching to a
    # numpy compressed zip archive file works.
    out = np.empty(nvals, dtype=list(zip(colnames, dtypes)))

    out["times"][:] = dat["pitch"].times[ok]
    out["pitch"][:] = dat["pitch"].vals[ok]
    out["roll"][:] = dat["roll"].vals[ok]
    out["alpang"][:] = dat["aoalpang"].vals[ok]
    out["betang"][:] = dat["aobetang"].vals[ok]
    out["alpha_sun"][:] = np.char.strip(dat["aoalpsun"].vals[ok]) == "SUN"
    out["beta_sun"][:] = np.char.strip(dat["aobetsun"].vals[ok]) == "SUN"
    out["kalman"][:] = (dat["aoacaseq"].vals[ok] == "KALM") & (
        dat["aopcadmd"].vals[ok] == "NPNT"
    )
    return out


get_fss_prim_data_cached = cache_file(get_fss_prim_data)


def add_pitch_roll_columns(dat, obc=True):
    """Add pitch, roll, pitch_err, roll_err columns to a table of FSS data

    Parameters
    ----------
    dat : numpy.ndarray or astropy.table.Table
        Table of FSS data with aoalpang and aobetang columns
    obc : bool
        Use OBC (True) or ground (False) FSS sun vector calibration.

    Returns
    -------
    dat : astropy.table.Table
        Table of FSS data with pitch_fss and roll_fss columns added.
    """
    dat = Table(dat)
    pitch_fss, roll_fss = get_pitch_roll_fss(dat["alpang"], dat["betang"], obc=obc)
    dat["pitch_fss"] = pitch_fss.astype(np.float32)
    dat["roll_fss"] = roll_fss.astype(np.float32)
    dat["roll_err"] = dat["roll_fss"] - dat["roll"]
    dat["pitch_err"] = dat["pitch_fss"] - dat["pitch"]

    return dat


def get_large_fss_roll_errors(dat, threshold=1.2):
    """Get FSS data with large roll errors.

    Returns a Table with columns::

        date: date of observation
        pitch: pitch (deg)
        roll: roll (deg)
        kalman: kalman mode (bool)
        roll_err: roll error (deg)

    Parameters
    ----------
    dat : astropy.table.Table
        Table of FSS data
    threshold : float
        Threshold for large roll error (deg)

    Returns
    -------
    out : astropy.table.Table
        Table of FSS data with large roll errors
    """
    droll = dat["roll_fss"] - dat["roll"]
    ok = (
        dat["alpha_sun"]
        & dat["beta_sun"]
        & (dat["pitch"] < 135)
        & (np.abs(droll) > threshold)
    )
    dok = dat[ok]
    dok["date"] = CxoTime(dok["times"]).date
    out = dok["date", "pitch", "roll", "kalman"]
    out["roll_err"] = dok["roll_fss"] - dok["roll"]
    out["pitch"].format = ".2f"
    out["roll"].format = ".2f"
    out["roll_err"].format = ".2f"
    return out


def get_fss_data_intervals(
    dat,
    mask=None,
    sun_presence=None,
    pitch_in_fov=None,
    pitch=None,
    dpitch=1.0,
    roll=None,
    droll=1.0,
    pitch_err=None,
    dpitch_err=1.0,
    roll_err=None,
    droll_err=1.0,
    dt_join=1000,
):
    """Get intervals of FSS data within boxes of pitch or roll

    Parameters
    ----------
    dat : numpy.ndarray or astropy.table.Table
        Table of FSS data
    mask : numpy.ndarray or None
        Boolean mask to apply to ``dat``.
    sun_presence : bool, None
        If not None select with matching sun presence.
    pitch_in_fov : bool, None
        If not None select data with pitch in the FSS FOV.
    pitch : float or None
        Select data within ``dpitch`` of ``pitch`` value if supplied (deg).
    dpitch : float
        Pitch width threshold (deg)
    roll : float or None
        Select data within ``droll`` of ``roll`` value if supplied (deg).
    droll : float
        Roll width threshold (deg)
    pitch_err : float or None
        Select data within ``dpitch_err`` of ``pitch_err`` value if supplied (deg).
    dpitch_err : float
        Pitch error outlier threshold (deg)
    roll_err : float or None
        Select data within ``droll_err`` of ``roll_err`` value if supplied (deg).
    droll_err : float
        Roll error outlier threshold (deg)
    dt_join : float
        Join intervals separated by less than ``dt_join`` (sec)

    Returns
    -------
    intervals : Table of (start, stop) tuples
        Intervals of FSS data that are outliers in pitch or roll
    """
    from cheta.utils import logical_intervals
    from kadi.events.models import fuzz_states

    ok = np.ones(len(dat), dtype=bool)
    if mask is not None:
        ok &= mask
    if sun_presence is not None:
        sp = dat["alpha_sun"] & dat["beta_sun"]
        ok &= sp if sun_presence else ~sp
    if pitch_in_fov is not None:
        pif = dat["pitch"] < 135
        ok &= pif if pitch_in_fov else ~pif
    if pitch is not None:
        ok &= np.abs(dat["pitch"] - pitch) < dpitch
    if roll is not None:
        ok &= np.abs(dat["roll"] - roll) < droll
    if pitch_err is not None:
        ok &= np.abs(dat["pitch_fss"] - dat["pitch"] - pitch_err) < dpitch_err
    if roll_err is not None:
        ok &= np.abs(dat["roll_fss"] - dat["roll"] - roll_err) < droll_err

    intervals = logical_intervals(dat["times"], ok, max_gap=33)
    intervals = fuzz_states(intervals, dt_join)
    intervals["duration"].format = ".1f"
    intervals["tstart"].format = ".1f"
    intervals["tstop"].format = ".1f"

    return intervals


def get_spm_pitch_roll():
    """Get pitch/roll values for Sun Position Monitor.

    From safe_mode_handbook for SPM::

        if ACA_Sx < KP.Region_Bound :
            a = KP.X_ellipse(1)
            b = KP.Y_ellipse(1)
        else:
            a = KP.X_ellipse(2)
            b = KP.Y_ellipse(2)

        (ACA_Sx/a)^2 + (ACA_Sy/b)^2 <= 1?

    From OBC PCAD code pcad_kons.ada::

        -- X-axis boundary of allowed region
            Region_Bound	: Float_Type :=
                0.0872; -- On-Orbit    SPR 571

        -- X-semi-axis of ellipse defining allowed region
            X_Ellipse		: Vector2_Type :=
                ( 1.0, 0.69466 ); -- On-Orbit     SPR 571

        -- Y-semi-axis of ellipse defining allowed region
            Y_Ellipse		: Vector2_Type :=
                ( 0.342, 0.2588 ); -- On-Orbit       <-- SPR 548
    """

    KP_Region_Bound = 0.0872
    KP_X_ellipse = [1.0, 0.69466]
    KP_Y_ellipse = [0.342, 0.2588]

    theta = np.linspace(0, 2 * np.pi, 10000)
    thetas = []
    pitchs = []
    rolls = []
    for i_constraint in 0, 1:
        a = KP_X_ellipse[i_constraint]
        b = KP_Y_ellipse[i_constraint]

        x = a * np.cos(theta)
        y = b * np.sin(theta)
        z = -np.sqrt((1 - x**2 - y**2).clip(0, 1))

        pitch = np.rad2deg(arccos_clip(x))
        roll = np.rad2deg(np.arctan2(-y, -z))

        ok = x < KP_Region_Bound if i_constraint == 0 else x >= KP_Region_Bound
        thetas.append(theta[ok])
        pitchs.append(pitch[ok])
        rolls.append(roll[ok])

    theta = np.hstack(thetas)
    pitch = np.hstack(pitchs)
    roll = np.hstack(rolls)
    # Sort by theta
    idx = np.argsort(theta)
    return pitch[idx], roll[idx]
