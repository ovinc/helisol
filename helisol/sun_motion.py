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


import numpy as np

from .general import CONSTANTS
from .general import _minus_pi_to_pi, _fraction_of_day, _day_time, _get_time
from .earth_motion import Earth



# ============================ Final calculations ============================


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
        """Init sun object from specific date/time.

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
        self.utc_time = _get_time(utc_time)
        self._calculate_sun_position()

    @property
    def right_ascension(self):
        """Right ascension in degrees"""
        ε = np.radians(self.earth.orbit.axial_tilt)
        λapp = np.radians(self.earth.apparent_longitude)
        α = np.arctan2(np.cos(ε) * np.sin(λapp), np.cos(λapp))
        return np.degrees(α)

    @property
    def declination(self):
        """Declination in degrees"""
        ε = np.radians(self.earth.orbit.axial_tilt)
        λapp = np.radians(self.earth.apparent_longitude)
        ẟ = np.arcsin(np.sin(λapp) * np.sin(ε))
        return np.degrees(ẟ)

    @property
    def equation_of_time(self):
        """Equation of time in degrees."""
        α = np.radians(self.right_ascension)
        γ0 = np.radians(self.earth.orbit.spring_longitude)
        m = np.radians(self.earth.average_motion)
        return np.degrees(_minus_pi_to_pi(α - γ0 - m))

    @property
    def hourly_angle(self):
        """Hourly angle in degrees"""
        H = 2 * np.pi * (_fraction_of_day(time) - 0.5) - equation_of_time + longitude
        return np.degrees(H)

    def _calculate_sun_position(self):
        """Calculate sun angles (azimuth, height, etc.) from time"""

        self.earth = Earth(utc_time=self.utc_time)

        lat, long = [x * np.pi / 180 for x in self.location]  # conversion to radians



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
