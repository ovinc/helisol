# About

**helisol** is a Python 3 package that calculates precisely the earth trajectory (tilt, longitude, true anomaly etc.) and apparent sun position (height above horizon, azimuth etc.) in the sky seen from any specified position on Earth and at any time. Sunset, solar noon, and sunrise times are also provided.

- `Earth` describes the motion of earth around the sun, irrespective of the location on the surface of earth

- `Sun` describes the apparent position of the sun in the sky at a given location on earth.

The package also includes tools to store and manipulate angles (`Angle` class) and times (`Time` class).

*SOON*:
- Sundial calculation and generation
- Improvements in precision (currently, order of 1-2 seconds in time)

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

a = Angle(degrees=30)
sin(a)
a.sin()  # equivalent to line above

b = Angle(radians=np.pi/4)
tan(b)

# Allowed operations
a + b
a - b
a * 2
2 * a
a / 2

# And combinations, e.g.
(a - b) / 4
cos(2 * a - b)

# access values in different units
a.radians
a.degrees
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
```


# Requirements

Python >= 3.6

*Packages*
- numpy


# Author

Olivier Vincent
(ovinc.py@gmail.com)

# Contributors

Gilbert Vincent (provided all equations)

# References

Meeus, J. *Calculs astronomiques à l’usage des amateurs*. (Société astronomique de France, 2014).


# License

GNU GPLv3, see *LICENSE* file
