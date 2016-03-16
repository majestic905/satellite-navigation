import math
import numpy as np
from translations import *

def minutes_since_spring_equinox(when):
    """
    :type when: datetime.datetime
    :param when: datetime.datetime object in UTC timezone

    :rtype: int
    :return: Number of minutes passed since spring equinox
    """
    b = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]
    return ((sum(b[0:when.month]) + when.day - 81) % 365)*1440 + (when.hour % 24)*60 + when.minute


def visible_from_ded_position(position, satellite):
    """
    :type position: numpy.array
    :param position: 1x3 numpy.array with position coordinates

    :type satellite: numpy.array
    :param satellite: 1x3 numpy.array with satellite coordinates

    :rtype: bool
    :return: True if the satellite is visible from current position
    """
    return position.dot(satellite - position) >= 0


def distance(position, satellite):
    return np.linalg.norm(position-satellite)


def distance_derivative(position, satellite, distance, spherical):
    if spherical:
        R, psi, lam = position
        x, y, z = satellite

        dq_dpsi = R/distance * (x*math.sin(psi)*math.cos(lam) + y*math.sin(psi)*math.sin(lam) - z*math.cos(psi))
        dq_dlam = R/distance * math.cos(psi)*(x*math.sin(lam) - y*math.cos(lam))

        return [dq_dpsi, dq_dlam]
    else:
        return (position-satellite) / distance
