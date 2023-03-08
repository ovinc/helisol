"""General functions and tools for the helisol package."""

# ----------------------------- License information --------------------------

# This file is part of the helisol python package.
# Copyright (C) 2023 Olivier Vincent

# The helisol package is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# The helisol package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with the helisol python package.
# If not, see <https://www.gnu.org/licenses/>


# ================================= Imports ==================================


import datetime
import numpy as np


# ================================ CONSTANTS =================================


CONSTANTS = {}

CONSTANTS['reference time'] = datetime.datetime(2000, 1, 1, 12, 0)

a = {}  # Coefficients for average motion
a['L'] = 280.46646, 36000.76983, 0.0003032   # J. MEEUS SAF
a['M'] = 357.52911, 35999.05029, -0.0001537  # (degrees)

CONSTANTS['average motion coefficients'] = a
CONSTANTS['nutation coefficients'] = (125.04, -1934.136)
CONSTANTS['aberration coefficients'] = -0.00569, -0.00478


# =============================== MISC. tools ================================


def _minus_pi_to_pi(x):
    """Angle x modulo 2*pi, between -pi and pi"""
    return np.arctan2(np.sin(x), np.cos(x))


def _fraction_of_day(time):
    """return 0 for midnight, 0.5 for noon"""
    return (time.hour * 3600 + time.minute * 60 + time.second) / (24 * 3600)


def _day_time(fraction_of_day):
    """Inverse function for fraction_of_day().

    Using hack from
    https://stackoverflow.com/questions/656297/python-time-timedelta-equivalent
    """
    Δt = datetime.timedelta(days=1) * fraction_of_day
    today = datetime.date.today()
    midnight = datetime.time()
    dtime = datetime.datetime.combine(today, midnight) + Δt
    return dtime.time()
