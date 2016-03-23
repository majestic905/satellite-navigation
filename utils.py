from math import sin, cos
import numpy as np
from earth import *
from translations import *

def minutes_since_spring_equinox(when):
    """
    :param when: datetime.datetime object in UTC timezone
    :return: Number of minutes passed since spring equinox
    """
    b = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]
    days = (sum(b[0:when.month]) + when.day - 81) % 365
    hours = when.hour % 24
    minutes = when.minute
    return days*1440 + hours*60 + minutes


def sat_visible(position, satellite):
    """
    :param position: [x, y, z] numpy.array with position rectangular coordinates
    :param satellite: [x, y, z] numpy.array with satellite rectangular coordinates
    :return: True if the satellite is visible from current position
    """
    return position.dot(satellite - position) >= 0


def sat_really_visible(position, satellite):
    """
    :param position: [x, y, z] numpy.array with position rectangular coordinates
    :param satellite: [x, y, z] numpy.array with satellite rectangular coordinates
    :return: True if the satellite is visible from current position
    """
    left = position.dot(satellite - position)/(np.linalg.norm(satellite) * np.linalg.norm(satellite - position))
    right = - np.sqrt(1 - Earth.R ** 2 / np.linalg.norm(position) ** 2)
    # TODO: check for accuracy
    return left > right


def distance(position, satellite):
    """
    :param position: [x, y, z] numpy.array with position rectangular coordinates
    :param satellite: [x, y, z] numpy.array with satellite rectangular coordinates
    :return: int distance between them in Euclidean metric
    """
    return np.linalg.norm(position-satellite)


def distance_derivative(position, satellite, distance, spherical):
    """
    :param position: [R, psi, lam] if spherical==True, [x, y, z] otherwise
    :param satellite: [x, y, z] numpy.array with satellite rectangular coordinates
    :param distance: int distance between them in Euclidean metric
    :param spherical: bool specifying the type of position coordinates
    :return: either [dq/dpsi, qd/dlam] numpy.array if spherical==True, or [dq/dx, dq/dy, dq/dz] otherwise
    """
    if spherical:
        R, psi, lam = position

        # TODO: figure out which one is correct. for now it seems working both ways O_o

        # x, y, z = satellite
        # dq_dpsi = R/distance * (x*sin(psi)*cos(lam) + y*sin(psi)*sin(lam) - z*cos(psi))
        # dq_dlam = R/distance * cos(psi)*(x*sin(lam) - y*cos(lam))

        x, y, z = spherical_to_rectangular(position)
        xs, ys, zs = satellite
        dq_dpsi = R/distance * ((xs-x)*sin(psi)*cos(lam) + (ys-y)*sin(psi)*sin(lam) - (zs-z)*cos(psi))
        dq_dlam = R/distance * cos(psi)*((xs-x)*sin(lam) - (ys-y)*cos(lam))

        return np.array([dq_dpsi, dq_dlam])
    else:
        return (position-satellite) / distance


def rho_dot(position, sat_position, speed, sat_speed, dist):
        return (position-sat_position).dot(speed-sat_speed) / dist


def rho_dot_derivative(position, sat_position, speed, sat_speed, dist):
    args = [position, sat_position, speed, sat_speed, dist]
    return (speed-sat_speed)/dist - rho_dot(*args)*(position-sat_position)/(dist**2)
