"""Manage command line parsing for the helisol module."""

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

import argparse


descr = "Calculate sun position, sunset, sunrise, and solar noon."

parser = argparse.ArgumentParser(
    description=descr, formatter_class=argparse.RawTextHelpFormatter
)

msg = "Input location"

# The nargs='?' is to have a positional argument with a default value
parser.add_argument("location", type=str, nargs="?", help=msg)
args = parser.parse_args()
