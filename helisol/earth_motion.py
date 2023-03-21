"""Calculate earth motion around sun"""

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

from .general import CONSTANTS
from .general import _minus_pi_to_pi, Time


class EarthOrbit:
    """Store general characteristics of earth orbit"""

    def __init__(self, utc_time=None):
        """Init earth orbit object from specific date/time.

        Parameters
        ----------
        utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        EarthOrbit('9am')
        EarthOrbit(utc_time='June 10 10:08:44')
        """
        self.time = Time(utc=utc_time)

    def __repr__(self):
        return f'Earth Orbit on {self.utc_time.date()} at {self.utc_time.time().strftime("%H:%M:%S")} (UTC)'

    @property
    def excentricity(self):
        """Excentricity of earth's orbit"""
        return 0.016709 - 0.000042 * self.time.julian_centuries

    @property
    def axial_tilt(self):
        """Earth's axial tilt (degrees)"""
        t = self.time.julian_centuries
        return (23 + 26 / 60 + (21.5 - 46.8 * t) / 3600)

    @property
    def spring_longitude(self):
        """Longitude between perigee and spring (degrees)"""
        a = CONSTANTS['average motion coefficients']
        t = self.time.julian_centuries
        return sum([(a['L'][i] - a['M'][i]) * t**i for i in range(3)])

    @property
    def nutation_with_aberration(self):
        """(delta_lambda). (degrees)"""
        t = self.time.julian_centuries
        nut0, nut1 = CONSTANTS['nutation coefficients']
        abr0, abr1 = CONSTANTS['aberration coefficients']
        nut = (nut0 + nut1 * t) * np.pi / 180
        return (abr0 + abr1 * np.sin(nut))


class Earth:

    def __init__(self, utc_time=None):
        """Init earth object from specific date/time.

        Parameters
        ----------
        utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        Earth('9am')
        Earth(utc_time='June 10 10:08:44')
        """
        self.time = Time(utc=utc_time)
        self.orbit = EarthOrbit(utc_time=utc_time)

    @property
    def average_motion(self):
        """Average motion / perigee. t time in centuries. In degrees"""
        a = CONSTANTS['average motion coefficients']
        t = self.time.julian_centuries
        m = sum([a['M'][i] * t**i for i in range(3)]) * np.pi / 180
        return np.degrees(_minus_pi_to_pi(m))

    def anomaly_radians(self, iteration=5):
        """Iterative way of calculating excentric anomaly (in radians)"""
        m = np.radians(self.average_motion)
        if iteration == 0:
            return m
        else:
            u = self.anomaly_radians(iteration=iteration - 1)
            e = self.orbit.excentricity
            return m + e * np.sin(u)

    @property
    def true_anomaly(self):
        """theta: Kepler angle (omega * t). in degrees """
        u = self.anomaly_radians(iteration=CONSTANTS['anomaly iterations'])
        e = self.orbit.excentricity
        nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(u / 2))
        return np.degrees(nu)

    @property
    def longitude(self):
        """Longitude (in degrees)"""
        λ0 = self.earth.orbit.spring_longitude
        λ = λ0 + self.true_anomaly
        return λ

    @property
    def apparent_longitude(self):
        """Longitude (in degrees)"""
        λ = self.longitude
        Δλ = self.earth.orbit.nutation_with_aberration
        return λ + Δλ


