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


from copy import copy
import numpy as np

from .general import CONSTANTS, Angle, Time, refraction
from .general import sin, cos, tan
from .earth_motion import Earth
from .locations import Location


class Sun:
    """Sun description irrespective of location on earth"""

    def __init__(self, utc_time=None):
        """Init sun object at specific date/time.

        Parameters
        ----------
        - utc_time: datetime or str (default None, i.e. current time)

        Examples
        --------
        Sun('9am')
        Sun(utc_time='June 10 10:08:44')
        """
        self.update(utc_time=utc_time)

    def __repr__(self):
        date = self.time.utc.date()
        time = self.time.utc.time().strftime("%H:%M:%S")
        a = f'Sun on {date} at {time} (UTC)\n'

        ẟ = self.declination.degrees
        alpha = self.right_ascension.degrees
        eqt = self.equation_of_time.degrees
        b = f'Declination [{ẟ:.1f}°], Right Asc. [{alpha:.1f}°], EQT [{eqt:.1f}°]'

        return a + b

    def update(self, utc_time=None):
        """Update at current time (default) or given UTC time."""
        self.time = Time(utc_time=utc_time)
        self.earth = Earth(utc_time=utc_time)

    @property
    def angular_diameter(self):
        d = 1_392_000
        l = self.earth.distance
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


class SunObservation:
    """Observation of sun at specific location on earth."""

    def __init__(self, location=None, utc_time=None):
        """Init sun observation object from specific location and date/time.

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
        self.sun = Sun(utc_time=utc_time)

    def __repr__(self):

        lat = self.latitude.degrees
        long = self.longitude.degrees
        date = self.time.utc.date()
        time = self.time.utc.time().strftime("%H:%M:%S")
        a = f'Sun observation from ({lat}°, {long}°) on {date} at {time} (UTC)'

        height = self.height.degrees
        app_h = self.apparent_height.degrees
        azimuth = self.azimuth.degrees
        b = f'\nAzimuth [{azimuth:.1f}°], Height [{height:.1f}°], App. Height [{app_h:.1f}°]'

        moments = self.sunrise.utc, self.noon.utc, self.sunset.utc
        sunrise, noon, sunset = [m.strftime("%H:%M:%S") for m in moments]
        c = f'\nSunrise {sunrise} Noon {noon} Sunset {sunset}'

        return a + b + c

    @property
    def hourly_angle(self):
        x = self.time.fraction_of_day
        a = Angle(radians=(2 * np.pi * (x - 0.5)))
        return a - self.sun.equation_of_time + self.longitude

    @property
    def height(self):
        """Height of sun above horizon."""
        ẟ = self.sun.declination
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
        ẟ = self.sun.declination
        H = self.hourly_angle
        lat = self.latitude
        a = cos(ẟ) * sin(H)
        b = cos(ẟ) * cos(H) * sin(lat) - sin(ẟ) * cos(lat)
        return Angle.arctan2(a, b)

    # --------------------- Sunrise, solar noon, sunset ----------------------

    def _noon(self):
        """Raw solar noon, calculated for current time."""
        eqt = self.sun.equation_of_time
        long = self.longitude
        noon_frac = (eqt.radians - long.radians) / (2 * np.pi) + 0.5
        return Time(utc_time=self.time.utc.date(), fraction_of_day=noon_frac)

    def _sunrise(self):
        """Sunrise"""
        ẟ = self.sun.declination
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
            old_obs = self._converge_to_event(event, iteration=iteration - 1)
            event_time = getattr(old_obs, raw_event_name)()
            return SunObservation(old_obs.location, utc_time=event_time)

    @property
    def noon(self):
        return self._converge_to_event(event='noon')._noon()

    @property
    def sunrise(self):
        return self._converge_to_event(event='sunrise')._sunrise()

    @property
    def sunset(self):
        return self._converge_to_event(event='sunset')._sunset()

    def _actual_event(self,
                      event='sunrise',
                      refract=True,
                      point='center',
                      obstacle=Angle(0),
                      precision=Angle(degrees=0.01),
                      print_details=False):
        """Find actual moment of sunset or sunrise (with refraction / obstacles).
        (iteratively)

        Parameters
        ----------
        - event: 'sunrise' or 'sunset'
        - refract: (default: True): take into account refraction or not
        - point: consider 'center', 'top', or 'bottom' of the sun
        - obstacle: angular height of obstacle masking the sun
        - precision: which (angular) tolerance to consider matching heights
        - print_details: print info on the iteration / convergence process
        """
        coeffs = {'sunrise': -1, 'sunset': 1}
        point_coeff = {'center': 0, 'bottom': -1, 'top': 1}
        c = coeffs[event]
        p = point_coeff[point]

        # Initial guess: ~7 minutes before sunrise (or after sunset) at horizon
        f0 = self.sunrise.fraction_of_day + c * 5e-3
        obs_search = copy(self)

        def match_heights(f):
            """Function that returns 0 sun point matches given height"""
            time = Time(obs_search.time, fraction_of_day=f)
            obs_search.update(utc_time=time)

            h0 = obs_search.height
            diam = obs_search.sun.angular_diameter

            h = h0 + p * diam / 2
            if refract:
                h += refraction(obs_search.height)

            return (h - obstacle).degrees

        def mvtime(f0, step=1e-4, max_it=1e4):
            """Iterative search."""
            tolerance = precision.degrees
            i = 0
            f = f0
            m0 = match_heights(f0)
            direction = np.sign(m0) * c

            stop = False
            found = False
            max_iterations = False

            while not (found or stop or max_iterations):
                m = match_heights(f)
                if abs(m) < tolerance:
                    found = True
                elif np.sign(m) == - np.sign(m0):
                    stop = True
                elif i > max_it:
                    max_iterations = True
                else:
                    i += 1
                    f += direction * step

            return {'value': f,
                    'found': found,
                    'stopped': stop,
                    'iterations': i,
                    'max iterations': max_iterations,
                    'residual': m}

        def manage_result(results):
            results['iterations'] = total_iterations
            if print_details:
                print(results)
            f = results['value']
            return Time(utc_time=self.time, fraction_of_day=f)

        f = f0
        total_iterations = 0

        for step in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
            results = mvtime(f, step=step)
            f = results['value']
            total_iterations += results['iterations']
            if results['found']:
                return manage_result(results)
        else:
            if results['found']:
                return manage_result(results)
            else:
                msg = f'Impossible to converge {event} search within tolerance'
                raise RuntimeError(msg)

    def actual_sunrise(self, *args, point='top', **kwargs):
        """Find actual moment of sunrise (with refraction / obstacles), iteratively.

        Parameters
        ----------
        - refract: (default: True): take into account refraction or not
        - point: consider 'center', 'top', or 'bottom' of the sun
        - obstacle: angular height of obstacle masking the sun
        - precision: which (angular) tolerance to consider matching heights
        - print_details: print info on the iteration / convergence process
        """
        return self._actual_event(event='sunrise', *args, point=point, **kwargs)

    def actual_sunset(self, *args, point='top', **kwargs):
        """Find actual moment of sunrise (with refraction / obstacles), iteratively.

        Parameters
        ----------
        - refract: (default: True): take into account refraction or not
        - point: consider 'center', 'top', or 'bottom' of the sun
        - obstacle: angular height of obstacle masking the sun
        - precision: which (angular) tolerance to consider matching heights
        - print_details: print info on the iteration / convergence process
        """
        return self._actual_event(event='sunset', *args, point=point, **kwargs)
