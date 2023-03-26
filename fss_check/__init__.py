# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Ska_helpers is a collection of utilities for the Ska3 runtime environment.
"""
from ska_helpers.version import get_version

__version__ = get_version(__package__)

from .check_fss import *  # noqa: F403, F401
