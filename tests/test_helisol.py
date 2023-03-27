"""Test helisol module with pytest"""

from pathlib import Path
import pandas as pd

import helisol
from helisol import Sun, Time, SunObservation, Angle, refraction

DATA_FOLDER = Path(helisol.__file__).parent.parent / 'data'


# ================================ Misc tools ================================

def predict_event(row, event='sunrise'):
    obs = SunObservation((47, 2), utc_time=row['Date'])
    if event in ('sunrise', 'sunset'):
        ppty = 'actual_' + event
        event_time = getattr(obs, ppty)(point='center')
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


def test_sun():
    """Test sun object"""
    sun = Sun('Aug 10, 2023, 12:00')
    assert round(sun.declination.degrees, 1) == 15.6
    assert round(sun.equation_of_time.degrees, 1) == 1.4
    assert round(sun.right_ascension.degrees, 1) == 140.1


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
    sunset = obs.actual_sunset(point='center').rounded_to('minute')
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


def test_ephemerides():
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
