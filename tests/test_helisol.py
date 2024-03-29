"""Test helisol module with pytest"""

from pathlib import Path
import pandas as pd
import numpy as np

import helisol
from helisol import Earth, Sun, Time, SunObservation, Distance, astronomical_unit
from helisol import Angle, AngleArray, AngleFromDegrees, AngleFromRadians
from helisol import AngleFromMinutes, AngleFromSeconds
from helisol import refraction

DATA_FOLDER = Path(helisol.__file__).parent.parent / 'data'


# ================================ Misc tools ================================


def predict_event(row, event='sunrise'):
    obs = SunObservation((47, 2), utc_time=row['Date'])
    if event in ('sunrise', 'sunset'):
        ppty = 'actual_' + event
        event_time = getattr(obs, ppty)(point='center', precision=0.002)
    else:
        event_time = getattr(obs, event)
    return event_time.utc.time()


def event_diff(row, event='sunrise'):
    column_eph = event.capitalize()
    column_th = column_eph + ' (predicted)'
    event_eph = Time(row[column_eph]).utc
    event_th = Time(row[column_th]).utc
    return (event_th - event_eph).total_seconds()


# ================================== Tests ===================================


def test_angle():
    """Test basic angle operation"""
    a = Angle(30)
    b = Angle(radians=(np.pi * 60 / 180))
    c = Angle(hms=(2, 0, 0))
    assert round((a + b - c).cos(), 3) == 0.5


def test_specialized_angle():
    """Test faster Angle implementations"""
    a = AngleFromDegrees(30)
    b = AngleFromRadians(np.pi / 6)
    c = AngleFromMinutes(1.5)
    d = 90 * AngleFromSeconds(1)
    assert round((a - b + c - d).seconds, 2) == 0


def test_angle_array():
    """Test AngleArray class."""
    a = AngleArray(degrees=[44.5, 2.5, 3], minutes=[30, 30, 0])
    assert round(a.tan()[0], 6) == 1


def test_distance():
    """Test Distance class"""
    d = Distance(au=1)
    assert round(d.km * 1000 / astronomical_unit, 6) == 1


def test_aphelion():
    """Test aphelion calculation from distance function"""
    aph_time = Earth(2023).orbit.aphelion(resolution='hour')
    assert aph_time.utc.month == 7
    assert aph_time.utc.day == 6


def test_perihelion():
    """Test perihelion calculation from distance function"""
    peri_time = Earth(2023).orbit.perihelion(resolution='hour')
    assert peri_time.utc.month == 1
    assert peri_time.utc.day == 4


def test_sun():
    """Test sun object"""
    sun = Sun('Aug 10, 2023, 12:00')
    assert round(sun.declination.degrees, 1) == 15.6
    assert round(sun.equation_of_time.degrees, 1) == 1.4
    assert round(sun.right_ascension.degrees, 1) == 140.1
    assert round(sun.angular_diameter.degrees, 2) == 0.53


def test_sun_observation():
    """Test calculation of sun height and azimuth"""
    obs = SunObservation(location=(42.4, -76.5))
    obs.update(utc_time='Jan 6, 2023, 4:25:03pm')
    assert round(obs.azimuth.degrees, 2) == -11.84
    assert round(obs.height.degrees, 2) == 24.24


def test_sunset():
    """Test calculation of sunset time"""
    obs = SunObservation(location=(42.4, -76.5), utc_time='2023-1-6, 16:25:03')
    assert obs.sunset.utc.hour == 21
    assert obs.sunset.utc.minute == 43


def test_actual_sunset():
    obs = SunObservation(location=(47, 2), utc_time='2023-03-15')
    sunset = obs.actual_sunset(point='center', precision=0.002).rounded_to('minute')
    assert sunset.utc.hour == 17
    assert sunset.utc.minute == 56


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


def test_ephemerides_march2023():
    """Check that predicted values coincide with Ephemerides"""
    file = DATA_FOLDER / 'Ephemerides_March2023_47N_2E.tsv'
    df = pd.read_csv(file, sep='\t')
    for event in ('sunrise', 'noon', 'sunset'):
        pred = lambda x: predict_event(x, event)
        diff = lambda x: event_diff(x, event)
        column_pred = event.capitalize() + ' (predicted)'
        column_diff = event.capitalize() + ' (diff)'
        df[column_pred] = df.apply(pred, axis=1)
        df[column_diff] = df.apply(diff, axis=1)
    assert (df['Noon (diff)'].abs() < 1).all()
    assert (df['Sunset (diff)'].abs() < 60).all()
    assert (df['Sunrise (diff)'].abs() < 60).all()


def test_ephemerides_aug2023():
    file = DATA_FOLDER / 'Ephemerides_Aug2023.tsv'
    df = pd.read_csv(file, sep='\t')

    # Right_ascension --------------------------------------------------------

    def analyze_asc(row):
        utc_time = row['Date'] + ' ' + row['Time']
        sun = Sun(utc_time=utc_time)
        h, m, s = [int(x) for x in row['Right Ascension'].split(':')]
        asc_eph = Angle(hms=(h, m, s))
        asc_th = sun.right_ascension
        return (asc_eph.degrees - asc_th.degrees) * 240  # in time-seconds

    df['Right Ascension (diff)'] = df.apply(analyze_asc, axis=1)
    max_dev_asc_time_s = df['Right Ascension (diff)'].abs().describe()['max']

    # Declination ------------------------------------------------------------

    def analyze_decl(row):
        utc_time = row['Date'] + ' ' + row['Time']
        sun = Sun(utc_time=utc_time)
        deg, mm, ss = [int(x) for x in row['Declination'].split(':')]
        decl_eph = Angle(degrees=deg, minutes=mm, seconds=ss)
        decl_th = sun.declination
        return (decl_eph.degrees - decl_th.degrees) * 3600  # in arc-seconds

    df['Declination (diff)'] = df.apply(analyze_decl, axis=1)
    max_dev_decl_arc_s = df['Declination (diff)'].abs().describe()['max']

    # Equation of time -------------------------------------------------------

    def analyze_eqt(row):
        utc_time = row['Date'] + ' ' + row['Time']
        sun = Sun(utc_time=utc_time)
        m, s = [int(x) for x in row['Equation of Time'].split(':')]
        eqt_eph = Angle(hms=(0, m, s))
        eqt_th = sun.equation_of_time
        return (eqt_eph.degrees - eqt_th.degrees) * 240

    df['Equation of Time (diff)'] = df.apply(analyze_eqt, axis=1)
    max_dev_eqt_time_s = df['Equation of Time (diff)'].abs().describe()['max']

    # Earth-sun distance -----------------------------------------------------

    def analyze_dist(row):
        utc_time = row['Date'] + ' ' + row['Time']
        sun = Sun(utc_time=utc_time)
        return sun.earth.distance.au

    df['Distance to Earth (predicted)'] = df.apply(analyze_dist, axis=1)
    df['Distance to Earth (diff, km)'] = (df['Distance to Earth'] - df['Distance to Earth (predicted)']) * astronomical_unit / 1000
    df.plot(y=['Distance to Earth', 'Distance to Earth (predicted)'])
    max_dev_dist_km = df['Distance to Earth (diff, km)'].describe()['max']

    # Final tests ------------------------------------------------------------

    assert max_dev_asc_time_s < 1     # Right asc. error < 1 seconds (in time)
    assert max_dev_decl_arc_s < 2.5   # Declination error < 2.5 arc-seconds
    assert max_dev_eqt_time_s < 1     # EQT error < 1 seconds (in time)
    assert max_dev_dist_km < 1600    # Distance earth-sun, error < 10000 km
