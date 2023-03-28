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
from dateutil.parser import parse
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

planet_coeffs = {'Venus 1': (134e-5, 351.52, 22518.443),
                 'Venus 2': (153e-5, 253.14, 45036.886),
                 'Jupiter': (200e-5, 157.23, 32964.467),
                 'Moon': (180e-5, 297.85, 445267.112),
                 'Long period': (196e-5, 252.08, 20.190)}

CONSTANTS['perturbations'] = planet_coeffs

astronomical_unit = 149_597_870_700  # in meters


# ============================== MISC. classes ===============================


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
        test_degs = []
        for x in (degrees, minutes, seconds):
            try:
                test_deg = bool(abs(x) > 0)
            except ValueError:  # array
                test_deg = (abs(x) > 0).any()
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
        a = "helisol.Angle\n"
        try:
            b = f"""{sign}{round_deg}°{minutes}'{seconds:.2f}" """
            c = f"[{sign}{h}h{m}m{s:.2f}s]\n"
        except TypeError:
            b = f"""{sign}{round_deg}°{minutes}'{np.round(seconds, 2)} """
            c = f"[{sign}{h}h{m}m{np.round(s, 2)}s]\n"
        d = f'{self.degrees} [°]\n'
        e = f'{self.radians} [rad]'
        return a + b + c + d + e

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
        return self._radians

    @radians.setter
    def radians(self, value):
        self._radians = value
        self._degrees = np.degrees(value)

    @property
    def degrees(self):
        return self._degrees

    @degrees.setter
    def degrees(self, value):
        self._degrees = value
        self._radians = np.radians(value)

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


# ===================== Convenience functions on angles ======================


def sin(angle):
    return angle.sin()


def cos(angle):
    return angle.cos()


def tan(angle):
    return angle.tan()


def cotan(angle):
    return angle.cotan()


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
        h = true_height.degrees
        y = Angle(degrees=h + (10.3 / (h + 5.11)))
        return Angle(minutes=1.02 / tan(y))
    else:
        h0 = apparent_height.degrees
        y0 = Angle(degrees=h0 + 7.31 / (h0 + 4.4))
        return Angle(minutes=cotan(y0))
