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

from .general import Time
from .general import _minus_pi_to_pi, _fraction_of_day, _day_time
from .earth_motion import Earth


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
        self.latitude_radians, self.longitude_radians = [np.radians(x) for x in self.location]

        self.time = Time(utc=utc_time)
        self.earth = Earth(utc_time=utc_time)

    def __repr__(self):
        lat_deg, long_deg = self.location
        a = f'Sun seen from ({lat_deg}°, {long_deg}°) on {self.time.utc.date()} at {self.time.utc.time().strftime("%H:%M:%S")} (UTC)'
        b = f'\nHeight {self.height:.2f}, Azimuth {self.azimuth:.2f}'

        sunrise, noon, sunset = [t.strftime("%H:%M:%S") for t in (self.rise, self.noon, self.set)]
        c = f'\nSunrise {sunrise} Noon {noon} Sunset {sunset}'
        return a + b + c

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
        time = self.time.utc.time()
        eqt = np.radians(self.equation_of_time)
        H = 2 * np.pi * (_fraction_of_day(time) - 0.5) - eqt + self.longitude_radians
        return np.degrees(H)

    @property
    def height(self):
        """Height of sun above horizon in radians"""
        ẟ = np.radians(self.declination)
        lat = self.latitude_radians
        H = np.radians(self.hourly_angle)
        h = np.arcsin(np.sin(ẟ) * np.sin(lat) + np.cos(ẟ) * np.cos(lat) * np.cos(H))
        return np.degrees(h)

    @property
    def azimuth(self):
        """Azimuth of sun with respect to South"""
        ẟ = np.radians(self.declination)
        H = np.radians(self.hourly_angle)
        lat = self.latitude_radians
        a = np.cos(ẟ) * np.sin(H)
        b = np.cos(ẟ) * np.cos(H) * np.sin(lat) - np.sin(ẟ) * np.cos(lat)
        A = np.arctan2(a, b)
        return np.degrees(A)

    @property
    def solar_noon(self):
        """Solar noon in fraction of day (0 for midnight, 0.5 for noon)"""
        eqt = np.radians(self.equation_of_time)
        long = self.longitude_radians
        return (eqt - long) / (2 * np.pi) + 0.5

    @property
    def sunrise(self):
        ẟ = np.radians(self.declination)
        lat = self.latitude_radians
        return self.solar_noon - np.arccos(-np.tan(ẟ) * np.tan(lat)) / (2 * np.pi)

    @property
    def sunset(self):
        return 2 * self.solar_noon - self.sunrise

    @property
    def noon(self):
        return _day_time(self.solar_noon)

    @property
    def rise(self):
        return _day_time(self.sunrise)

    @property
    def set(self):
        return _day_time(self.sunset)
