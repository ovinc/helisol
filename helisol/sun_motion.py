"""Calculate sun motion."""

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

from .general import CONSTANTS
from .general import _minus_pi_to_pi, _fraction_of_day, _day_time


# ================= Functions depending on time in centuries =================


def _excentricity(t):
    """Excentricity of earth's orbit; t time in centuries"""
    return 0.016709 - 0.000042 * t


def _obliquity(t):
    """Earth's axial tilt; t time in centuries"""
    return (23 + 26 / 60 + (21.5 - 46.8 * t) / 3600) * np.pi / 180


def _longitude_of_perigee(t):
    """Longitude perigee / gamma. t time in centuries. output in radians."""
    a = CONSTANTS['average motion coefficients']
    return sum([(a['L'][i] - a['M'][i])  * t**i for i in range(3)]) * np.pi / 180


def _average_sun_motion(t):
    """Average motion / perigee. t time in centuries. output in radians."""
    a = CONSTANTS['average motion coefficients']
    m = sum([a['M'][i] * t**i for i in range(3)]) * np.pi / 180
    return _minus_pi_to_pi(m)


def _nutation_with_aberration(t):
    """(delta_lambda). t: time in centuries"""
    nut0, nut1 = CONSTANTS['nutation coefficients']
    abr0, abr1 = CONSTANTS['aberration coefficients']
    nut = (nut0 + nut1 * t) * np.pi / 180
    return (abr0 + abr1 * np.sin(nut)) * np.pi / 180


# ========= Functions that depend on average motion and excentricity =========


def _anomaly(average_motion, excentricity, iteration=5):
    """Iterative way of calculating excentric anomaly"""
    if iteration == 0:
        return average_motion
    else:
        u = _anomaly(average_motion, excentricity=excentricity, iteration=iteration - 1)
        return average_motion + excentricity * np.sin(u)


def _true_anomaly(average_motion, excentricity):
    """theta: Kepler angle (omega * t). """
    e = excentricity
    u = _anomaly(average_motion=average_motion, excentricity=e)
    return 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(u / 2))


# ======== Functions that depend on apparent longitude and obliquity =========


def _right_ascension(apparent_longitude, obliquity):
    """Right ascension in radians"""
    return np.arctan2(np.cos(obliquity) * np.sin(apparent_longitude),
                      np.cos(apparent_longitude))


def _declination(apparent_longitude, obliquity):
    """Declination in radians"""
    return np.arcsin(np.sin(apparent_longitude) * np.sin(obliquity))


# ============================ Final calculations ============================


def _equation_of_time(average_motion, right_ascension, longitude_of_perigee):
    """Equation of time, in radians"""
    return _minus_pi_to_pi(right_ascension - longitude_of_perigee - average_motion)


def _hourly_angle(time, equation_of_time, longitude):
    """Hourly angle in radians"""
    return 2 * np.pi * (_fraction_of_day(time) - 0.5) - equation_of_time + longitude


def _height_above_horizon(declination, hourly_angle, latitude):
    """Height of sun above horizon in radians"""
    return np.arcsin(np.sin(declination) * np.sin(latitude) + np.cos(declination) * np.cos(latitude) * np.cos(hourly_angle))


def _azimuth_south(declination, hourly_angle, latitude):
    """Azimuth of sun with respect to South"""
    a = np.cos(declination) * np.sin(hourly_angle)
    b = np.cos(declination) * np.cos(hourly_angle) * np.sin(latitude) - np.sin(declination) * np.cos(latitude)
    return np.arctan2(a, b)


# ================================ MAIN CLASS ================================


class Sun:

    def __init__(self, location, utc_time=None):
        """Init sun object from spcific date/time.

        Parameters
        ----------
        location: tuple (or iterable) (latitude, longitude) in degrees
        utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        loc = (42.4, -76.5)
        Sun(loc, '9am')
        Sun(location=loc, utc_time='June 10 10:08:44')
        """
        self.location = location
        self.update(utc_time=utc_time)

    def __repr__(self):
        lat_deg, long_deg = self.location
        a = f'Sun seen from ({lat_deg}°, {long_deg}°) on {self.utc_time.date()} at {self.utc_time.time().strftime("%H:%M:%S")} (UTC)'
        b = f'\nHeight {self.height:.2f}, Azimuth {self.azimuth:.2f}'

        sunrise, noon, sunset = [t.strftime("%H:%M:%S") for t in (self.rise, self.noon, self.set)]
        c = f'\nSunrise {sunrise} Noon {noon} Sunset {sunset}'
        return a + b + c

    @staticmethod
    def _get_time(utc_time=None):
        """Return a datetime object from user input"""
        if utc_time is None:
            return datetime.datetime.utcnow()
        else:
            return parse(str(utc_time), yearfirst=True)  # str is in case a datetime object is passed

    def update(self, utc_time=None):
        """Update sun position at current time or time set by user.

        Parameters
        ----------
        location: tuple (or iterable) (latitude, longitude) in degrees
        utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        sun.update()
        sun.update('17:56:05')
        sun.update(utc_time=datetime.datetime(2022, 10, 11, 15, 04))
        """
        self.utc_time = self._get_time(utc_time)
        self._calculate_sun_position()

    def _calculate_sun_position(self):
        """Calculate sun angles (azimuth, height, etc.) from time"""
        lat, long = [x * np.pi / 180 for x in self.location]  # conversion to radians

        elapsed_time = self.utc_time - CONSTANTS['reference time']
        time_days = elapsed_time.total_seconds() / (24 * 3600)
        time_centuries = time_days / (365.25 * 100)

        e = _excentricity(t=time_centuries)
        ε = _obliquity(t=time_centuries)
        m = _average_sun_motion(t=time_centuries)
        λ0 = _longitude_of_perigee(t=time_centuries)
        Δλ = _nutation_with_aberration(t=time_centuries)

        λ = λ0 + _true_anomaly(average_motion=m, excentricity=e)
        λapp = λ + Δλ

        right_asc = _right_ascension(apparent_longitude=λapp, obliquity=ε)
        ẟ = _declination(apparent_longitude=λapp, obliquity=ε)

        eqt = _equation_of_time(average_motion=m,
                                right_ascension=right_asc,
                                longitude_of_perigee=λ0)

        H = _hourly_angle(time=self.utc_time.time(),
                          equation_of_time=eqt,
                          longitude=long)

        height = _height_above_horizon(declination=ẟ, hourly_angle=H, latitude=lat)
        azimuth = _azimuth_south(declination=ẟ, hourly_angle=H, latitude=lat)

        self.height = height * 180 / np.pi
        self.azimuth = azimuth * 180 / np.pi

        solar_noon = (eqt - long) / (2 * np.pi) + 0.5
        sunrise = solar_noon - np.arccos(-np.tan(ẟ) * np.tan(lat)) / (2 * np.pi)
        sunset = 2 * solar_noon - sunrise

        self.noon = _day_time(solar_noon)
        self.rise = _day_time(sunrise)
        self.set = _day_time(sunset)
