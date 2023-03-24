# About

**helisol** is a Python 3 package that calculates precisely the earth trajectory (tilt, longitude, true anomaly etc.) and apparent sun position (height above horizon, azimuth etc.) in the sky seen from any specified position on Earth and at any time. Sunset, solar noon, and sunrise times are also provided.

- `Earth` describes the motion of earth around the sun, irrespective of the location on the surface of earth

- `Sun` describes the apparent position of the sun in the sky at a given location on earth.

The package also includes tools to store and manipulate angles (`Angle` class) and times (`Time` class), and a function to calculate refraction effects (`refraction()`).

*SOON*:
- Sundial calculation and generation
- Improvements in precision (currently, order of a few seconds in time)

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

## Sun viewed from given location on earth

```python
from helisol import Sun

# Get current position of the sun
# NOTE: location can be a tuple of coords, a location name (if stored in the
# JSON database), or a location object (see below).
sun = Sun(location=(42.4, -76.5))
print(sun)  # Some info (azimuth height, sunrise etc. is printed here)

# Update position to current time
sun.update()
print(sun)

# Update position to specified time and date
sun.update(utc_time='Jan 6, 2023, 4:25:03pm')
print(sun.height)   # Height above horizon, in degrees
print(sun.azimuth)  # azimuth with respect to south in degrees

# It is possible to specify time upon instantiation directly:
sun = Sun(location=(42.4, -76.5), utc_time='2023-1-6, 16:25:03')
print(sun.sunrise, sun.noon, sun.sunset)
```

## Angles

```python
from helisol import Angle, sin, cos, tan

a = Angle(degrees=30.7)
sin(a)
a.sin()  # equivalent to line above

b = Angle(radians=np.pi/4)
tan(b)

c = Angle(degrees=3, minutes=45, seconds=10)
cos(c)

d = Angle.arctan(1)  # pi/4

# Allowed operations
a + b
a - b
a * 2
2 * a
a / 2
- a

# And combinations, e.g.
(a - b) / 4
cos(2 * a - b)

# access values in different units
a.radians
a.degrees

# Modify value in-place
a.degrees = 60
a.radians = np.pi / 4

# Read-only info on minutes and seconds of angle
a.minutes
a.seconds
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
From helisol import Location, Sun

# Load existing location and use it to instantiate a Sun object
location = Location.load('Home')
sun = Sun(location)
# equivalently:
sun = Sun('Home')

# It is possible to configure a default location in config.py (default 'Home')
# so that one can do just
Sun()

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
