import numpy as np
from math import sin, cos, atan2, sqrt


# TODO: CHECK THE FUCKING ACCURACY. WHERE THE MINUS IS SUPPOSED TO BE?
def greenwich_to_equatorial(greenwich, gamma):
    matrix = np.array([[ cos(gamma),  -sin(gamma), 0],
                       [sin(gamma),  cos(gamma), 0],
                       [         0,           0, 1]])
    return matrix.dot(greenwich)


# TODO: CHECK THE FUCKING ACCURACY. WHERE THE MINUS IS SUPPOSED TO BE?
def equatorial_to_greenwich(equatorial, gamma):
    matrix = np.array([[cos(gamma), sin(gamma), 0],
                       [-sin(gamma),  cos(gamma), 0],
                       [          0,          0, 1]])
    return matrix.dot(equatorial)


def spherical_to_rectangular(spherical):
    """
    :param spherical: numpy.array([R, psi, lam])
    :return: numpy.array([x, y, z])
    """
    R, psi, lam = spherical
    return R*np.array([cos(psi)*cos(lam),
                       cos(psi)*sin(lam),
                       sin(psi)])


def rectangular_to_spherical(rectangular):
    """
    :param rectangular: numpy.array([x, y, z])
    :return: numpy.array([R, psi, lam])
    """
    x, y, z = rectangular
    R = np.linalg.norm([x, y, z])
    psi = atan2(z, sqrt(x**2 + y**2))
    lam = atan2(y, x)
    return np.array([R, psi, lam])
