import numba as nb
import numpy as np

@nb.jit(nopython=True, cache=True)
def wgs84_to_ecef(lla):
    coords_ecef = np.zeros(3)
    
    a = 6378137.0
    b = 6356752.314245179
    ba = 0.9966471893352525

    psi = np.arctan( np.tan(lla[0]) * ba )
    r = a * np.cos(psi) + lla[2] * np.cos(lla[0])

    coords_ecef[0] = r * np.cos(lla[1])
    coords_ecef[1] = r * np.sin(lla[1])
    coords_ecef[2] = b * np.sin(psi) + lla[2] * np.sin(lla[0])

    return coords_ecef


@nb.njit(cache=True)
def TEC(tau: float, delta: float):
    """Transformation from :ref:`sec:F-C` to the :ref:`sec:F-E`

    :param tau: Longitude (radians) from the Greenwich meridian (tau is positive if the vehicle position is east of the Greenwich meridian)
    :type tau: float
    :param delta: Latitude (radians) from the equator (delta is positive if the vehicle location is on the northern hemisphere)
    :type delta: float
    :return: TEC
    :rtype: numpy_array"""

    sang = np.sin(np.array([tau, delta]))
    cang = np.cos(np.array([tau, delta]))

    return np.array([[-sang[1]*cang[0], -sang[1]*sang[0], cang[1]],
                     [-sang[0], cang[0], 0.],
                     [-cang[1]*cang[0], -cang[1]*sang[0], -sang[1]]])
