# About

**helisol** is a Python 3 package that calculates precisely the apparent sun position (height above horizon, azimuth) in the sky seen from any specified position on Earth and at any time. Sunset, solar noon, and sunrise times are also provided.

*SOON*: Sundial calculation and generation.

# Install

```bash
pip install helisol
```

# Usage

```python
from helisol import Sun

# Get current position of the sun
sun = Sun(location=(42.4, -76.5))
print(sun)  # All info (azimuth height, sunrise etc. is printed here)

# Update position to current time
sun.update()
print(sun)

# Update position to specified time and date
sun.update(utc_time='Jan 6, 2023, 4:25:03pm')
print(sun.height)   # Height above horizon, in degrees
print(sun.azimuth)  # azimuth with respect to south in degrees

# It is possible to specify time upon instantiation directly:
sun = Sun(location=(42.4, -76.5), utc_time='2023-1-6, 16:25:03')
print(sun.rise, sun.noon, sun.set)
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
