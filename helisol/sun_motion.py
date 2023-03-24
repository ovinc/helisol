"""Calculate sun position"""

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


import numpy as np

from .general import CONSTANTS, Angle, Time, refraction
from .general import sin, cos, tan
from .earth_motion import Earth
from .locations import Location


class Sun:

    def __init__(self, location=None, utc_time=None):
        """Init sun object from specific date/time.

        Parameters
        ----------
        - location: Location name from JSON database,
                    [or] custom Location object
                    [or] iterable (latitude, longitude) in degrees
                    if None (default), use default location

        - utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        loc = (42.4, -76.5)
        Sun(loc, '9am')
        Sun(location=loc, utc_time='June 10 10:08:44')
        """
        self.location = Location.parse(location)
        self.latitude, self.longitude = [Angle(degrees=x) for x in self.location.coords]
        self.update(utc_time=utc_time)

    def update(self, utc_time=None):
        """Update at current time (default) or given UTC time."""
        self.time = Time(utc_time=utc_time)
        self.earth = Earth(utc_time=utc_time)

    def __repr__(self):
        lat = self.latitude.degrees
        long = self.longitude.degrees
        date = self.time.utc.date()
        time = self.time.utc.time().strftime("%H:%M:%S")
        a = f'Sun seen from ({lat}°, {long}°) on {date} at {time} (UTC)'

        height = self.height.degrees
        azimuth = self.azimuth.degrees
        b = f'\nHeight {height:.2f}, Azimuth {azimuth:.2f}'

        moments = self.sunrise.utc, self.noon.utc, self.sunset.utc
        sunrise, noon, sunset = [m.strftime("%H:%M:%S") for m in moments]
        c = f'\nSunrise {sunrise} Noon {noon} Sunset {sunset}'
        return a + b + c

    @property
    def earth_distance(self):
        l0 = 149_500_000
        e = self.earth.orbit.excentricity
        nu = self.earth.true_anomaly
        return l0 * (1 - e**2) / (1 + e * cos(nu))

    @property
    def angular_diameter(self):
        d =  1_392_000
        l = self.earth_distance
        return Angle.arctan(d / l)

    @property
    def right_ascension(self):
        ε = self.earth.orbit.axial_tilt
        λapp = self.earth.apparent_longitude
        return Angle.arctan2(cos(ε) * sin(λapp), cos(λapp))

    @property
    def declination(self):
        ε = self.earth.orbit.axial_tilt
        λapp = self.earth.apparent_longitude
        return Angle.arcsin(sin(λapp) * sin(ε))

    @property
    def equation_of_time(self):
        alpha = self.right_ascension
        gamma0 = self.earth.orbit.spring_longitude
        m = self.earth.average_motion
        Δaberr = self.earth.orbit.correc_aberration
        Δψ = self.earth.orbit.correc_nutation
        ε = self.earth.orbit.axial_tilt

        eqt = (alpha - gamma0 - m - Δaberr - Δψ * cos(ε))
        eqt.minus_pi_to_pi()

        return eqt

    @property
    def hourly_angle(self):
        x = self.time.fraction_of_day
        a = Angle(radians=(2 * np.pi * (x - 0.5)))
        return a - self.equation_of_time + self.longitude

    @property
    def height(self):
        """Height of sun above horizon."""
        ẟ = self.declination
        lat = self.latitude
        H = self.hourly_angle
        return Angle.arcsin(sin(ẟ) * sin(lat) + cos(ẟ) * cos(lat) * cos(H))

    @property
    def apparent_height(self):
        """Apparent height including refraction"""
        h = self.height
        return h + refraction(true_height=h)

    @property
    def azimuth(self):
        """Azimuth of sun with respect to South"""
        ẟ = self.declination
        H = self.hourly_angle
        lat = self.latitude
        a = cos(ẟ) * sin(H)
        b = cos(ẟ) * cos(H) * sin(lat) - sin(ẟ) * cos(lat)
        return Angle.arctan2(a, b)

    # --------------------- Sunrise, solar noon, sunset ----------------------

    def _noon(self):
        """Raw solar noon, calculated for current time."""
        eqt = self.equation_of_time
        long = self.longitude
        noon_frac = (eqt.radians - long.radians) / (2 * np.pi) + 0.5
        return Time(utc_time=self.time.utc.date(), fraction_of_day=noon_frac)

    def _sunrise(self):
        """Sunrise"""
        ẟ = self.declination
        lat = self.latitude
        noon_frac = self._noon().fraction_of_day
        rise_frac = noon_frac - np.arccos(-tan(ẟ) * tan(lat)) / (2 * np.pi)
        return Time(utc_time=self.time.utc.date(), fraction_of_day=rise_frac)

    def _sunset(self):
        """Sunset"""
        set_frac = 2 * self._noon().fraction_of_day - self._sunrise().fraction_of_day
        return Time(utc_time=self.time.utc.date(), fraction_of_day=set_frac)

    def _converge_to_event(self, event='noon', iteration=CONSTANTS['sunset iterations']):
        """Iterative way of conver ging sun towards sunrise, noon, or sunset"""
        raw_event_name = '_' + event
        if iteration == 0:
            return self
        else:
            old_sun = self._converge_to_event(event, iteration=iteration - 1)
            event_time = getattr(old_sun, raw_event_name)()
            return Sun(old_sun.location, utc_time=event_time)

    @property
    def noon(self):
        return self._converge_to_event(event='noon')._noon()

    @property
    def sunrise(self):
        return self._converge_to_event(event='sunrise')._sunrise()

    @property
    def sunset(self):
        return self._converge_to_event(event='sunset')._sunset()
