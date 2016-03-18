import math
import numpy as np
from earth import *
from translations import *

def minutes_since_spring_equinox(when):
    """
    :param when: datetime.datetime object in UTC timezone
    :return: Number of minutes passed since spring equinox
    """
    b = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]
    return ((sum(b[0:when.month]) + when.day - 81) % 365)*1440 + (when.hour % 24)*60 + when.minute


def visible_from_ded_position(position, satellite):
    """
    :param position: 1x3 numpy.array with position eq coordinates
    :param satellite: 1x3 numpy.array with satellite eq coordinates
    :return: True if the satellite is visible from current position
    """
    return position.dot(satellite - position) >= 0


def really_visible_from_ded_position(position, satellite):
    """
    :param position: 1x3 numpy.array with position eq coordinates
    :param satellite: 1x3 numpy.array with satellite eq coordinates
    :return: True if the satellite is visible from current position
    """
    left = position.dot(satellite - position)/(np.linalg.norm(satellite) * np.linalg.norm(satellite - position))
    right = - np.sqrt(1 - Earth.R ** 2 / np.linalg.norm(position) ** 2)
    # TODO: check for accuracy
    # print(left, right, left > right)
    return left > right


def distance(position, satellite):
    """
    :param position: 1x3 numpy.array with position eq coordinates
    :param satellite: 1x3 numpy.array with satellite eq coordinates
    :return: distance between them in Euclidean metric
    """
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
