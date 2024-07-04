import numpy as np

from datetime import datetime, timedelta


def generate_dates(start_date_str, end_date_str):
    """

    generates list of dates between given start date and end date including both.
 
    Parameters
    ----------
 
    start_date_str : str
        fmt: "YYYYMMDD"

    end_date_str : str
        fmt: "YYYYMMDD"
  
    Returns
    -------

    out: list

    """

    date_format = "%Y%m%d"  # Adjust the format to match your date strings
    start_date = datetime.strptime(start_date_str, date_format)
    end_date = datetime.strptime(end_date_str, date_format)

    date_list = []
    current_date = start_date

    while current_date <= end_date:
        date_list.append(current_date.strftime("%Y%m%d"))
        current_date += timedelta(days=1)

    return date_list


def spaceAngle(theta1, phi1, theta2, phi2, rad=False):
    """
      calculates angular distance between (theta1, phi) & (theta2, phi2).
      all inputs should be in the units deg. if rad is set FALSE.
      outputs are in the units deg, if rad is set FALSE.
 
      Parameters
      ----------
  
      theta1 : float or array
	  zenith - 1	

      phi1 : float or array
          azimuth - 1

      theta2: float or array
          zenith - 2

      phi2: float or array
          azimuth - 2 
      Returns
      -------
      out: array 
          angular distance 
    """

    if rad == False:
        theta1 = np.radians(theta1)
        phi1 = np.radians(phi1)

        theta2 = np.radians(theta2)
        phi2 = np.radians(phi2)

        return np.degrees(
            np.arccos(
                np.cos(theta1) * np.cos(theta2) * np.cos(phi1 - phi2)
                + np.sin(theta1) * np.sin(theta2)
            )
        )

    else:
        return np.arccos(
            np.cos(theta1) * np.cos(theta2) * np.cos(phi1 - phi2)
            + np.sin(theta1) * np.sin(theta2)
        )


def distance_based_cut(CoreX, CoreY, log10NKGSize):

    distFromMuSt = np.sqrt((CoreX - -61.91) ** 2 + (CoreY - -19.35) ** 2)
    NKGsizeEdges = np.array(
        [0, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 20]
    )
    distValues = np.array(
        [0, 12, 24, 34, 44, 52, 62, 74, 86, 98, 112, 126, 126, 126, 0]
    )

    distcut = distValues[np.digitize(log10NKGSize, NKGsizeEdges) - 1]
    return distFromMuSt < distcut
