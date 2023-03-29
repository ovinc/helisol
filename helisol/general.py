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
import math
from dateutil.parser import parse
from functools import total_ordering
import numpy as np


# ================================ CONSTANTS =================================


CONSTANTS = {}

CONSTANTS['reference time'] = datetime.datetime(2000, 1, 1, 12, 0)

a = {}  # Coefficients for average motion
a['L'] = 280.46646, 36000.76983, 0.0003032   # J. MEEUS SAF
a['M'] = 357.52911, 35999.05029, -0.0001537  # (degrees)

CONSTANTS['average motion coefficients'] = a
CONSTANTS['nutation coefficients'] = (125.04, -1934.136)

CONSTANTS['anomaly iterations'] = 3
CONSTANTS['sunset iterations'] = 2

planet_coeffs = {'Venus 1': (351.52, 22518.443),
                 'Venus 2': (253.14, 45036.886),
                 'Jupiter': (157.23, 32964.467),
                 'Moon': (297.85, 445267.112),
                 'Long period': (252.08, 20.190),
                 'H': (353.40, 65928.7155)}

# amplitudes for longitude (left) and sun-earth radius (right)
planet_amplitudes = {'Venus 1': (134e-5, 5.43e-6),
                     'Venus 2': (153e-5, 15.75e-6),
                     'Jupiter': (200e-5, 16.27e-6),
                     'Moon': (180e-5, 30.76e-6),
                     'Long period': (196e-5, 0),
                     'H': (0, 9.27e-6)}

CONSTANTS['planet perturbation coefficients'] = planet_coeffs
CONSTANTS['planet perturbation amplitudes'] = planet_amplitudes

astronomical_unit = 149_597_870_700  # in meters


# ============================== MISC. classes ===============================


@total_ordering
class Angle:
    """Store angles and retrieve them in degrees or radians.

    Note: can also store arrays of angles.
    """
    def __init__(self,
                 degrees=0, minutes=0, seconds=0,
                 hms=None,
                 radians=None):
        """Input angle in degrees (+ minutes/seconds), or in radians, or
        in time-like format (h, m, s).

        Notes:
        - By default, value is considered to be decimal degrees, e.g. Angle(10.7)
        - degrees, hours, minutes, seconds can be floats and are added
          (but it's preferrable to use ints when using  deg, min, sec instead
           of decimal degrees)
        - same for h, m, s
        - ValueError raised if input has mixed units between °, hms, rad
        """
        # ATTENTION the nature of the tests below can significantly reduce
        # speed of some operations, e.g. calculating sunsets, sunrises etc.
        test_degs = []
        for x in (degrees, minutes, seconds):
            try:
                test_deg = bool(abs(x) > 0)
            except ValueError:  # if array is passed, consider it as an input
                test_deg = True
            finally:
                test_degs.append(test_deg)

        test_input_deg = any(test_degs)
        test_input_hms = (hms is not None)
        test_input_rad = (radians is not None)
        if sum([test_input_deg, test_input_hms, test_input_rad]) > 1:
            raise ValueError('Cannot specify angle in more than one set of units')

        if test_input_rad:
            # note: degrees etc. are set automatically thanks to setter
            self.radians = radians
        elif test_input_hms:
            h, m, s = hms
            self.degrees = 15 * (h + m / 60 + s / 3600)
        else:
            self.degrees = degrees + minutes / 60 + seconds / 3600

    def __repr__(self):
        if (np.array(self.degrees) < 0).all():
            sign = '-'
        elif (np.array(self.degrees) < 0).any():
            sign = '+/-'
        else:
            sign = ''
        round_deg = np.abs(np.int_(self.degrees))
        h, m, s = [abs(x) for x in self.hms]
        minutes = np.abs(self.minutes)
        seconds = np.abs(self.seconds)
        a = "helisol.Angle"
        try:  # angle with single value
            b = f"""\n{sign}{round_deg}°{minutes}'{seconds:.2f}" """
            c = f"[{sign}{h}h{m}m{s:.2f}s]\n"
        except TypeError:  # array
            b = " (array of angles)"
            c = "\n"
        d = f'{self.degrees} [°]\n'
        e = f'{self.radians} [rad]'
        return a + b + c + d + e

    def __eq__(self, other):
        return (self.degrees == other.degrees)

    def __lt__(self, other):
        return (self.degrees < other.degrees)

    def __add__(self, other):
        return Angle(degrees=self.degrees + other.degrees)

    def __sub__(self, other):
        return Angle(degrees=self.degrees - other.degrees)

    def __neg__(self):
        return Angle(degrees=-self.degrees)

    def __mul__(self, other):
        return Angle(degrees=self.degrees * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Angle(degrees=self.degrees / other)

    @property
    def radians(self):
        return self._degrees * np.pi / 180  # much faster than np.radians() for single values

    @radians.setter
    def radians(self, value):
        self._degrees = value * 180 / np.pi  # much faster than np.degrees() for single values

    @property
    def degrees(self):
        return self._degrees

    @degrees.setter
    def degrees(self, value):
        self._degrees = value

    @property
    def hms(self):
        """Only for info, not settable"""
        h_decimal = self.degrees / 15
        s = np.sign(h_decimal)
        h_abs = np.abs(h_decimal)
        h = np.int_(h_decimal)
        m = np.int_(s * 60 * np.remainder(h_abs, 1))
        s = s * 60 * np.remainder(h_abs * 60, 1)
        return (h, m, s)

    @property
    def minutes(self):
        """Only for info, not settable"""
        s = np.sign(self.degrees)
        deg_abs = np.abs(self.degrees)
        # if s in outside of int_(), returns floats instead of ints
        return np.int_(s * 60 * np.remainder(deg_abs, 1))

    @property
    def seconds(self):
        """Only for info, not settable"""
        s = np.sign(self.degrees)
        deg_abs = np.abs(self.degrees)
        return s * 60 * np.remainder(deg_abs * 60, 1)

    def sin(self):
        return np.sin(self.radians)

    def cos(self):
        return np.cos(self.radians)

    def tan(self):
        return np.tan(self.radians)

    def cotan(self):
        return 1 / np.tan(self.radians)

    def minus_pi_to_pi(self):
        """Angle modulo 2*pi, between -pi and pi"""
        self.radians = np.arctan2(self.sin(), self.cos())

    @classmethod
    def arcsin(cls, value):
        return Angle(radians=np.arcsin(value))

    @classmethod
    def arccos(cls, value):
        return Angle(radians=np.arccos(value))

    @classmethod
    def arctan(cls, value):
        return Angle(radians=np.arctan(value))

    @classmethod
    def arctan2(cls, x1, x2):
        return Angle(radians=np.arctan2(x1, x2))


class AngleFromDegrees(Angle):
    """Shortcut for faster execution of Angle when known to instantiate from degrees"""

    def __init__(self, value=0):
        self.degrees = value


class AngleFromMinutes(Angle):
    """Shortcut for faster execution of Angle when known to instantiate from arcminutes"""

    def __init__(self, value=0):
        self.degrees = value / 60


class AngleFromSeconds(Angle):
    """Shortcut for faster execution of Angle when known to instantiate from arcseconds"""

    def __init__(self, value=0):
        self.degrees = value / 3600


class AngleFromRadians(Angle):
    """Shortcut for faster execution of Angle when known to instantiate from radians"""

    def __init__(self, value=0):
        self.radians = value


# ===================== Convenience functions on angles ======================


def sin(angle):
    return angle.sin()


def cos(angle):
    return angle.cos()


def tan(angle):
    return angle.tan()


def cotan(angle):
    return angle.cotan()


# =========================== Distance management ============================


@total_ordering
class Distance:

    def __init__(self, m=0, km=None, au=None):
        """Init distance with meters, kilometers or astronomical units.

        Data can be single values or arrays.
        """
        try:
            test_input_m = bool(abs(m) > 0)
        except ValueError:  # if array is passed, consider it as an input
            test_input_m = True
        test_input_km = (km is not None)
        test_input_au = (au is not None)
        if sum([test_input_m, test_input_km, test_input_au]) > 1:
            raise ValueError('Cannot specify angle in more than one set of units')

        if test_input_km:
            self.km = km
        elif test_input_au:
            self.au = au
        else:
            self.m = m

    def __repr__(self):
        return f'helisol.Distance\n{self.m}[m]\n{self.km}[km]\n{self.au}[A.U.]'

    def __eq__(self, other):
        return (self.m == other.m)

    def __lt__(self, other):
        return (self.m < other.m)

    def __add__(self, other):
        return Distance(m=self.m + other.m)

    def __sub__(self, other):
        return Distance(m=self.m - other.m)

    def __neg__(self):
        return Distance(m=-self.m)

    def __mul__(self, other):
        return Distance(m=self.m * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Distance(m=self.m / other)

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, value):
        self._m = value
        self._km = value / 1000
        self._au = value / astronomical_unit

    @property
    def km(self):
        return self._km

    @km.setter
    def km(self, value):
        self._km = value
        self._m = value * 1000
        self._au = self._m / astronomical_unit

    @property
    def au(self):
        return self._au

    @au.setter
    def au(self, value):
        self._au = value
        self._m = value * astronomical_unit
        self._km = self._m / 1000


# =========================== Date/Time management ===========================


class Time:
    """Store time in various useful formats"""

    def __init__(self, utc_time=None, fraction_of_day=None):
        """Init time object.

        Parameters
        ----------
        - utc_time: datetime or str (default None, i.e. current time)

        - fraction_of_day: if specified, overrides time from the given fraction
                           of day (0.5 for Noon) [keeps input date]
        """
        self.utc = self._parse_time(utc_time)
        self._update()

        if fraction_of_day is not None:
            self.fraction_of_day = fraction_of_day

    def __repr__(self):
        return str(self.utc)

    def _update(self):
        """Recalculate quantities if UTC time has changed."""
        self.elapsed = self.utc - CONSTANTS['reference time']
        self.days = self.elapsed.total_seconds() / (24 * 3600)
        self.julian_years = self.days / 365.25
        self.julian_centuries = self.julian_years / 100

    def _parse_time(self, utc_time=None):
        """Return a datetime object from user input"""
        if utc_time is None:
            return datetime.datetime.utcnow()
        else:
            return parse(str(utc_time), yearfirst=True)  # str is in case a datetime object is passed

    @property
    def fraction_of_day(self):
        """return 0 for midnight, 0.5 for noon"""
        time = self.utc.time()
        return (time.hour * 3600 + time.minute * 60 + time.second) / (24 * 3600)

    @fraction_of_day.setter
    def fraction_of_day(self, value):
        """Inverse function for fraction_of_day().

        Kind of using hack from
        https://stackoverflow.com/questions/656297/python-time-timedelta-equivalent
        """
        Δt = datetime.timedelta(days=1) * value
        date = self.utc.date()
        midnight = datetime.time()
        self.utc = datetime.datetime.combine(date, midnight) + Δt
        self._update()

    def rounded_to(self, unit='second'):
        """Return rounded version of self.

        Parameters
        ----------
        - unit (str): 'second' or 'minute'
        """
        if unit == 'second':
            old_time = self.utc
            new_time = old_time.replace(microsecond=0)
            if old_time.microsecond >= 5e5:
                new_time += datetime.timedelta(seconds=1)

        elif unit == 'minute':
            old_time = self.rounded_to('second').utc
            new_time = old_time.replace(second=0)
            if old_time.second >= 30:
                new_time += datetime.timedelta(minutes=1)

        return Time(utc_time=new_time)


# ================================ Refraction ================================


def refraction_saemundsson_math(true_height):
    """Refraction angle from true height (Saemundsson 1986).

    Uses math.tan (faster for single values)

    Parameters
    ----------
    true_height: angle in degrees

    Output
    ------
    Refraction angle in arc-minutes
    """
    y = true_height + (10.3 / (true_height + 5.11))
    return 1.02 / math.tan(y * np.pi / 180)


def refraction_saemundsson_np(true_height):
    """Refraction angle from true height (Saemundsson 1986)

    Uses np.tan (for arrays)

    Parameters
    ----------
    true_height: angle in degrees

    Output
    ------
    Refraction angle in arc-minutes
    """
    y = true_height + (10.3 / (true_height + 5.11))  # in degrees
    return 1.02 / np.tan(y * np.pi / 180)


def refraction_bennett_np(apparent_height=None):
    """Refraction angle from apparent height (Bennett 1982)

    Parameters
    ----------
    apparent_height: angle in degrees

    Output
    ------
    Refraction angle in arc-minutes
    """
    y0 = apparent_height + 7.31 / (apparent_height + 4.4)  # in degrees
    return 1 / np.tan(y0 * np.pi / 180)


def refraction(true_height=None, apparent_height=None):
    """Refraction angle, either from true height or apparent height.

    Parameters
    ----------
    - true_height: helisol.Angle
    - apparent_height: helisol.Angle

    If both true_height and apparent_height are provided, apparent_height
    is ignored.

    True height formula is from Saemundsson 1986
    Apparent height formula is from Bennett 1982

    Normally valid for 1010 bars of pressure and 10°C of temperature.

    Examples
    --------
    refraction(true_height=Angle(0))
    refraction(apparent_height=Angle(0))
    refraction(Angle(23))                          # true height by default
    refraction(true_height=Angle(degrees=23))
    refraction(apparent_height=Angle(minutes=66))
    """
    if true_height is not None:
        return AngleFromMinutes(refraction_saemundsson_np(true_height.degrees))
    else:
        return AngleFromMinutes(refraction_bennett_np(apparent_height.degrees))
