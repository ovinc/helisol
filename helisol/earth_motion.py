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
from .general import Angle, Time
from .general import sin, cos, tan


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
        self.time = Time(utc_time=utc_time)

    def __repr__(self):
        year = self.time.utc.date().year
        return f'Earth Orbit in {year}'

    @property
    def excentricity(self):
        """Excentricity of earth's orbit"""
        return 0.016709 - 0.000042 * self.time.julian_centuries

    @property
    def axial_tilt(self):
        """Earth's axial tilt (epsilon)"""
        t = self.time.julian_centuries
        eps0 = Angle(degrees=(23 + 26 / 60 + (21.5 - 46.8 * t) / 3600))
        d_eps = Angle(seconds=9.2 * cos(self.nutation))  # tilt nutation
        return eps0 + d_eps

    @property
    def spring_longitude(self):
        """Longitude between perigee and spring"""
        a = CONSTANTS['average motion coefficients']
        t = self.time.julian_centuries
        gamma0 = sum([(a['L'][i] - a['M'][i]) * t**i for i in range(3)])
        return Angle(degrees=gamma0)

    @property
    def nutation(self):
        """Ω"""
        t = self.time.julian_centuries
        nut0, nut1 = CONSTANTS['nutation coefficients']
        return Angle(degrees=(nut0 + nut1 * t))

    @property
    def correc_nutation(self):
        """Δψ"""
        Ω = self.nutation
        return Angle(seconds=-17.2 * sin(Ω))

    @property
    def correc_aberration(self):
        """Δaberr"""
        return Angle(seconds=-20.5)

    @property
    def correc_planets(self):
        """Perturbations from planets, moon, etc."""
        t = self.time.julian_centuries
        perturb = Angle()
        funcs = cos, cos, cos, sin, sin
        for planet_coeffs, func in zip(CONSTANTS['perturbations'].values(), funcs):
            ampl, a0, a1 = planet_coeffs
            ag = Angle(degrees=(a0 + a1 * t))
            corr = Angle(degrees=(ampl * func(ag)))
            perturb += corr
        return perturb


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
        self.update(utc_time=utc_time)

    def update(self, utc_time=None):
        self.time = Time(utc_time=utc_time)
        self.orbit = EarthOrbit(utc_time=utc_time)

    def __repr__(self):
        date = self.time.utc.date()
        time = self.time.utc.time().strftime("%H:%M:%S")
        return f'Earth on {date} at {time} (UTC)'

    @property
    def average_motion(self):
        """Average motion with respect to perigee"""
        a = CONSTANTS['average motion coefficients']
        t = self.time.julian_centuries
        m_deg = sum([a['M'][i] * t**i for i in range(3)])
        m = Angle(degrees=m_deg)
        m.minus_pi_to_pi()
        return m

    def anomaly(self, iteration=5):
        """Iterative way of calculating excentric anomaly"""
        m0 = self.average_motion
        if iteration == 0:
            return m0
        else:
            u = self.anomaly(iteration=iteration - 1)
            e = self.orbit.excentricity
            m_rad = m0.radians + e * sin(u)
            return Angle(radians=m_rad)

    @property
    def true_anomaly(self):
        """True anomaly (nu)"""
        u = self.anomaly(iteration=CONSTANTS['anomaly iterations'])
        e = self.orbit.excentricity
        return 2 * Angle.arctan(np.sqrt((1 + e) / (1 - e)) * tan(u / 2))

    @property
    def longitude(self):
        """Longitude (λ) with respect to spring"""
        λ0 = self.orbit.spring_longitude
        λs = λ0 + self.true_anomaly
        return λs

    @property
    def apparent_longitude(self):
        """Apparent longitude (including correction from nutation, aberration, planets etc.)"""
        λs = self.longitude
        Δψ = self.orbit.correc_nutation
        Δaberr = self.orbit.correc_aberration
        Δplanets = self.orbit.correc_planets
        return λs + Δψ + Δaberr + Δplanets
