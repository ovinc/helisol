"""Test helisol module with pytest"""

from helisol import Sun, Angle, refraction


def test_sun_position():
    """Test calculation of sun height and azimuth"""
    sun = Sun(location=(42.4, -76.5))
    sun.update(utc_time='Jan 6, 2023, 4:25:03pm')
    assert round(sun.azimuth.degrees, 2) == -11.84
    assert round(sun.height.degrees, 2) == 24.24


def test_sunset():
    """Test calculation of sunset time"""
    sun = Sun(location=(42.4, -76.5), utc_time='2023-1-6, 16:25:03')
    assert sun.sunset.utc.hour == 21
    assert sun.sunset.utc.minute == 43


def test_refraction_true_height():
    """Test calculation of refraction from true height"""
    r1 = refraction(true_height=Angle(degrees=0))
    r2 = refraction(Angle(minutes=-34))
    assert r1.minutes == 28
    assert r2.minutes == 34


def test_refraction_apparent_height():
    """Test calculation of refraction from true height"""
    r1 = refraction(apparent_height=Angle(degrees=0))
    r2 = refraction(apparent_height=Angle(minutes=34))
    assert r1.minutes == 34
    assert r2.minutes == 28
