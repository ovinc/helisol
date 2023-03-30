# About

**helisol** is a Python 3 package that calculates precisely the earth trajectory (tilt, longitude, true anomaly etc.) and absolute and apparent sun position (declination, right ascension, height above horizon, azimuth etc.) in the sky seen from any specified position on Earth and at any time. Sunset, solar noon, and sunrise times are also provided. It uses Meeus' equations including perturbations from planets on the earth orbit.

- `Earth` describes the motion of earth around the sun

- `Sun` describes the position of the sun seen from earth, irrespective of exact location on earth (declination, equation of time, etc.)

- `SunObservation` describes the position of the sun in the sky as observed at a given location on earth (azimuth, height, etc.).

The package also includes tools to store and manipulate angles (`Angle` class) and times (`Time` class), and a function to calculate refraction effects (`refraction()`).

See *ExampleShadow.ipynb* for an example on how to use the package to calculate the shadow of an object at different times of the year.

*SOON*:
- Sundial calculation and generation

**Precision** can be evaluated by comparison to published ephemerids. We have used ephemerids for March 2023 for sunset / sunrise and solar noon, and ephemerids for August 2023 for right ascension, declination and equation of time (see test file).

- Sunset, sunrise and solar noon are within the displayed precision of each other (1 second for noon, 1 minute for sunset/sunrise)

- Right ascension and equation of time are within the displayed precision of each other (1 second in time)

- Declination deviates less than 2.5 arc-seconds (0.0007°, or 0.17 seconds in time).

- Sun to Earth distance deviates less than 1600 km.


# Install

```bash
pip install helisol
```

# Usage

## Earth position and orbit

```python
from helisol import Earth

# Get current axial tilt and true anomaly of the earth
earth = Earth()
earth.true_anomaly
earth.orbit.axial_tilt

# Get properties at a given date
earth = Earth('Sept. 9, 2003')

# Update to current time
earth.update()

# Update to specific time
earth.update(utc_time='2023-03-21, 12:00')
```


## Sun absolute position

```python
from helisol import Sun

# Get current sun position
sun = Sun()
sun.declination
sun.right_ascension
sun.equation_of_time

# Access earth object associated with sun
sun.earth

# Get sun at specific UTC time and then update to current time
sun = Sun('March 3, 13:30')
sun.update()
```


## Sun viewed from given location on earth

```python
from helisol import SunObservation

# Get current position of the sun
# NOTE: location can be a tuple of coords, a location name (if stored in the
# JSON database), or a location object (see below).
obs = Sun(location=(42.4, -76.5))
print(obs)  # Some info (azimuth height, sunrise etc. is printed here)

# Update position to current time
obs.update()
print(obs)

# Update position to specified time and date
obs.update(utc_time='Jan 6, 2023, 4:25:03pm')
print(obs.height)   # Height above horizon, in degrees
print(obs.azimuth)  # azimuth with respect to south in degrees

# It is possible to specify time upon instantiation directly:
obs = Sun(location=(42.4, -76.5), utc_time='2023-1-6, 16:25:03')

# sunrise, noon (meridian), sunset [center of sun, no refraction]
print(obs.sunrise, obs.noon, obs.sunset)

# To include refraction and other options to calculate sunrise / sunset:
obs.actual_sunrise()  # [top of sun, with refraction]
obs.actual_sunrise(point='center')  # [center of sun, with refraction]
obs.actual_sunrise(refract=False, point='bottom')  # [sun bottom, no refract]
obs.actual_sunrise(obstacle=Angle(26))  # with an obstacle of 26° of height
# The method also exists for sunset:
obs.actual_sunset()
# (NOTE: options for the precision of the calculation also exists, see docstring)

# Access sun object associated with the observation:
obs.sun
```

## Angles

```python
from helisol import Angle, sin, cos, tan

a = Angle(degrees=30.7)
sin(a)
a.sin()  # equivalent to line above

b = Angle(radians=np.pi/4)
tan(b)

c = Angle(degrees=3, minutes=45, seconds=10)  # °, ', "
cos(c)

d = Angle(hms=(8, 44, 43))  # in hours, min, sec

e = Angle.arctan(1)  # pi/4

# Allowed operations
a + b
a - b
a * 2
2 * a
a / 2
- a
a < b
a == b
# and other comparisons

# And combinations, e.g.
(a - b) / 4
cos(2 * a - b)

# access values in different units (decimal)
a.radians
a.degrees

# Read-only info on minutes and seconds of angle
a.minutes
a.seconds

# Read_only info on time in hours (min/sec):
a.hms
```

If having performance issues with `Angle`, it is possible to increase speed by shortcutting tests of units during instantiation, using more specialized `Angle` subclasses:

```python
from helisol import AngleFromDegrees, AngleFromRadians
from helisol import AngleFromMinutes, AngleFromSeconds
a = AngleFromDegrees(30)
b = AngleFromRadians(np.pi / 6)
c = AngleFromMinutes(1.5)
d = 90 * AngleFromSeconds(1)
a - b + c - d
```

It is also possible to use arrays with the `AngleArray` class, e.g.
```python
from helisol import AngleArray
a = AngleArray(degrees=[44.5, 2.5, 3], minutes=[30, 30, 0])
a.tan()
```

## Distances

Distances are managed by the `Distance` class

```python
from helisol import Distance
d = Distance(au=2.5)  # 2.5 astronomical units
d.m
d.km
d.au
```


## Date / times


```python
from helisol import Time

time = Time(utc_time='2022, July 3, 6:12')
time.utc  # datetime with no timezone (implicitly UTC)
time.julian_years  # Julian years since Jan 1 2000, 12:00
time.fraction_of_day  # 0 for midnight, 0.5 for noon

# It is possible to manually adjust the fraction of day to set the time:
time.fraction_of_day = 0.5

# It is also possible to instantiate a Time object with a fraction of day
time = Time(utc_time='2023, March 10', fraction_of_day=0.25)  # 6am

# To obtain a rounded version of time:
time.rounded_to('second')
time.rounded_to('minute')
```

## Refraction

```python
from helisol import Angle, refraction

# When true height is at horizon level (29' approx.)
refraction(true_height=Angle(0))

# When apparent height is at horizon level (34' approx)
refraction(apparent_height=Angle(0))

# At arbitrary true height or apparent heights
refraction(Angle(23))                          # true height of 23 degrees
refraction(true_height=Angle(degrees=23))      # same thing
refraction(apparent_height=Angle(minutes=66))  # apparent height of 10'
```


## Generate tables of data

It is possible to generate tables (pandas DataFrames) containing sun position data, or sunset/sunrise data between selected dates and with selected frequency, with the functions:
- `generate_table()`
- `extend_table()`
- `sunset_table()`

(see docstrings for help and examples).


## Location management

It is possible to save/load location information with the `Location` class.

```python
From helisol import Location, SunObservation

# Load existing location and use it to instantiate a Sun object
location = Location.load('Home')
obs = SunObservation(location)
# equivalently:
obs = SunObservation('Home')

# It is possible to configure a default location in config.py (default 'Home')
# so that one can do just
SunObservation()

# Define custom location and save it in the database
# NOTE: it is possible to define elevation, although not used in helisol at
# the moment.
location = Location(name='Work', coords=(40.78, -73.97))
location.save()          # save in a non-shared file (excluded from version )
location.save('global')  # save in globals.json file, version controlled.

# Remove location from database:
location.remove('Work')
```


# Requirements

Python >= 3.7

*Packages*
- numpy
- pandas
- oclock


# Author

Olivier Vincent
(ovinc.py@gmail.com)

# Contributors

Gilbert Vincent (equations and advice)

# References

- Meeus, J., *Calculs astronomiques à l’usage des amateurs*, Société astronomique de France (2014).

- Saemundsson, T., *Astronomical Refraction*, Sky and Telescope, **72**, 70 (1986).

- Bennett, G. G., *The calculation of astronomical refraction in marine navigation*, Journal of Navigation, **35**, 255-259 (1982).


# License

GNU GPLv3, see *LICENSE* file
