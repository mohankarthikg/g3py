from astropy.coordinates import EarthLocation, AltAz, ICRS, Angle, SkyCoord
from astropy.time import Time
import astropy.units as u
import numpy as np


def ZenPhiToRADec(theta, phi, EvDate, EvTime1, EvTime2):
    """
    Converts from GRAPES-3 Local coordindates to International Celestial Reference System (ICRS).

    Parameters
    ----------

    theta : float /  array_like 
        units: deg.
    phi : float / array_like
        units: deg.
    EvDate : int / str / array_like
        fmt: YYYYMMDD
    EvTime1 : int / str / array_like 
        fmt: HHMMSS
    EvTime2 : int / str /  array_like
        fmt: 9 digits (GRAPES-3 DAQ resolution is ns)
    
    Returns
    -------
    RightAscension: float / array
        units : deg.
    Declination: float / array 
        units : deg.
    """

    grapes3loc = EarthLocation(lon=76.7 * u.deg, lat=11.4 * u.deg, height=2200 * u.m)

    theta = np.array(theta)
    phi = np.array(phi)

    EvTime1 = np.array(EvTime1)
    EvTime2 = np.array(EvTime2)
    EvDate = np.array(EvDate)

    EvTime2 = np.char.zfill(EvTime2.astype("str"), width=9)
    EvTime1 = np.char.zfill(EvTime1.astype("str"), width=6)
    EvDate = np.char.zfill(EvDate.astype("str"), width=8)

    iteration = zip(list(EvDate), list(EvTime1), list(EvTime2))

    datetime = np.array(
        [
            f"{EvD[:4]}-{EvD[4:6]}-{EvD[6:]}\
 {EvT1[:2]}:{EvT1[2:4]}:{EvT1[4:6]}.{EvT2}"
            for EvD, EvT1, EvT2 in iteration
        ]
    )

    obs_time = Time(datetime, format="iso", location=grapes3loc)

    coord = SkyCoord(
        alt=(90 - theta) * u.deg,
        az=((180 - phi) % 360) * u.deg,
        frame="altaz",
        obstime=obs_time,
        location=grapes3loc,
    )

    coord = coord.transform_to("icrs")

    return coord.ra.degree, coord.dec.degree


def RADecToAltAz(ra, dec, EvDate, EvTime1, EvTime2):
    """
    Converts from International Celestial Reference System (ICRS) to GRAPES-3 local coordinates.

    Parameters
    ----------

    ra : float or array_like
        RightAscension.
        units: deg

    dec : float or array_like
        Declination.
        units: deg

    EvDate : int or str or array_like
        fmt: "YYYYMMDD"

    EvTime1 : int or str or array_like
        fmt: "HHMMSS"

    EvTime2 : int or str or array_like
        fmt: 9 digits (GRAPES-3 DAQ resolution is ns)

    Returns
    -------

    azimuth : float or array of same dim. as input
        units: deg

    altitude :  float or array of same dim. as input
        units: deg   
    """

    grapes3loc = EarthLocation(lon=76.7 * u.deg, lat=11.4 * u.deg, height=2200 * u.m)

    ra = np.array(ra)
    dec = np.array(dec)

    EvTime1 = np.array(EvTime1)
    EvTime2 = np.array(EvTime2)
    EvDate = np.array(EvDate)

    EvTime2 = np.char.zfill(EvTime2.astype("str"), width=9)
    EvTime1 = np.char.zfill(EvTime1.astype("str"), width=6)
    EvDate = np.char.zfill(EvDate.astype("str"), width=8)

    iteration = zip(list(EvDate), list(EvTime1), list(EvTime2))

    datetime = np.array(
        [
            f"{EvD[:4]}-{EvD[4:6]}-{EvD[6:]}\
 {EvT1[:2]}:{EvT1[2:4]}:{EvT1[4:6]}.{EvT2}"
            for EvD, EvT1, EvT2 in iteration
        ]
    )

    obs_time = Time(datetime, format="iso", location=grapes3loc)

    coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        frame="icrs",
        obstime=obs_time,
        location=grapes3loc,
    )

    coord = coord.transform_to("altaz")

    return (180 - coord.az.degree) % 360, coord.alt.degree


