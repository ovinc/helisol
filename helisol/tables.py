"""Generate data tables."""

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


import datetime
import pandas as pd
from oclock import parse_time

from .general import Time
from .sun_motion import Sun
from .locations import Location


column_names = {
    'azimuth': 'Azimuth (°)',
    'height': 'Height (°)',
    'apparent_height': 'Apparent Height (°)',
    'declination': 'Declination (°)',
    'right_ascension': 'Right Ascension (°)',
    'equation_of_time': 'Equation of Time (°)',
}


def _generate_times(start, end, interval):
    """Generate list of helisol.Time objects between specified times."""
    t1 = Time(start)
    t2 = Time(end)
    dt = parse_time(interval)

    times = []
    t = t1.utc
    while t <= t2.utc:
        times.append(Time(t))
        t += dt

    return times


def _round_to_second(time):
    """Round time to nearest second."""
    if time.microsecond < 5e5:
        return time.replace(microsecond=0)
    return time.replace(microsecond=0) + datetime.timedelta(seconds=1)


def generate_table(location, start, end, interval, columns=column_names):
    """Generate table with solar data at regular intervals between two dates.

    Parameters
    ----------
     - location: Location name from JSON database,
                    [or] custom Location object
                    [or] iterable (latitude, longitude) in degrees
                    if None (default), use default location

    - start: datetime, helisol.Time object, str or other datetime information (UTC)

    - end: datetime, helisol.Time object, str or other datetime information (UTC)

    - interval: str in the format 'hh:mm:ss', including '::10' for 10 secs

    - columns (optional): specify which properties and columns will be included
                          (dict with sun attribute name as key and name of
                           corresponding column as value)

    Output
    ------
    Pandas DataFrame

    Examples
    --------
    generate_table(location=(47, 2),
                   start='March 1 12:00',
                   end='March 31 12:00',
                   interval='24::',
                   columns={'azimuth': 'Azimuth (°)', 'height': 'Height (°)'})
    """
    location = Location.parse(location)
    times = _generate_times(start=start, end=end, interval=interval)

    data = {}
    data['Date'] = [time.utc.date() for time in times]
    data['Time (UTC)'] = [time.utc.time() for time in times]

    suns = [Sun(location=location, utc_time=time) for time in times]

    for ppty, name in columns.items():
        data[name] = [getattr(sun, ppty).degrees for sun in suns]

    return pd.DataFrame(data)


def sunset_table(location, start, end):
    """Create able with sunrise, noon and sunsets between specified dates.

    Parameters
    ----------
     - location: Location name from JSON database,
                    [or] custom Location object
                    [or] iterable (latitude, longitude) in degrees
                    if None (default), use default location

    - start: datetime, helisol.Time object, str or other datetime information (UTC)

    - end: datetime, helisol.Time object, str or other datetime information (UTC)

    Output
    ------
    Pandas DataFrame
    """
    times = _generate_times(start=start, end=end, interval='24::')

    data = {}
    data['Date'] = [time.utc.date() for time in times]
    data['Time (UTC)'] = [time.utc.time() for time in times]

    suns = [Sun(location=location, utc_time=time) for time in times]

    for ppty in 'sunrise', 'noon', 'sunset':
        name = ppty.capitalize()
        data[name] = [_round_to_second(getattr(sun, ppty).utc).time() for sun in suns]

    return pd.DataFrame(data)


def extend_table(data, location, date_column='Date', time_column='Time (UTC)', columns=column_names):
    """Takes an existing table with date / time columns and adds columns with
    the corresponding solar data.

    Parameters
    ----------
    - data: pandas DataFrame

     - location: Location name from JSON database,
                    [or] custom Location object
                    [or] iterable (latitude, longitude) in degrees
                    if None (default), use default location

    - date_column: name of the column containing dates in the pandas dataframe

    - time_column: name of the column containing time in the pandas dataframe

    - columns (optional): specify which properties and columns will be included
                          (dict with sun attribute name as key and name of
                           corresponding column as value)

    Example
    -------
    data = pd.read_excel('Data.xls')
    extend_table(data, location=(47, 2))
    """
    def _calculate_angle(ppty, row):
        utc_time = datetime.datetime.combine(row[date_column], row[time_column])
        sun = Sun(location=location, utc_time=utc_time)
        return getattr(sun, ppty).degrees

    for ppty, name in columns.items():
        calc = lambda row: _calculate_angle(ppty, row)
        data[name] = data.apply(calc, axis=1)
