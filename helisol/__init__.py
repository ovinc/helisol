"""Helisol package init."""

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


from .general import Time, Angle, sin, cos, tan, cotan, refraction
from .earth_motion import Earth
from .sun_motion import Sun, SunObservation
from .tables import generate_table, sunset_table, extend_table
from .locations import Location

# from importlib.metadata import version (only for python 3.8+)
from importlib_metadata import version


__version__ = version('helisol')
__author__ = 'Olivier Vincent'
__license__ = 'GNU GPLv3'
