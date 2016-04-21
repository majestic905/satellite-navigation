import numpy as np
from math import sin, cos, radians

from utils import *
from translations import *
from earth import *


class Satellite:
    # eps, omega = 0, 0

    def __init__(self, omega_big, tau):
        self.omega_big = radians(omega_big)
        self.tau = tau
        self.i, self.R, self.w = 0, 0, 0

    def coordinates_equatorial(self, minutes):
        # u = self.omega + self.w*(t-self.tau)  # omega = 0
        u = self.w*(minutes - self.tau)

        in_orbit_plane = self.R*np.array([cos(u), sin(u), 0])

        around_x = np.array([[1,           0,            0],
                             [0, cos(self.i), -sin(self.i)],
                             [0, sin(self.i),  cos(self.i)]])

        around_z = np.array([[cos(self.omega_big), -sin(self.omega_big), 0],
                             [sin(self.omega_big),  cos(self.omega_big), 0],
                             [                  0,                    0, 1]])

        return around_z.dot(around_x).dot(in_orbit_plane)

    def coordinates_greenwich(self, minutes):
        gamma = Earth.w_sun*minutes + Earth.w_self*minutes
        # return greenwich_to_equatorial(self.coordinates_equatorial(minutes), gamma)
        return equatorial_to_greenwich(self.coordinates_equatorial(minutes), gamma)

    def speed_equatorial(self, minutes):
        u = self.w*(minutes - self.tau)

        return self.R * self.w * \
               np.array([-sin(u)*cos(self.omega_big) - cos(u)*sin(self.omega_big)*cos(self.i),
                         -sin(u)*sin(self.omega_big) + cos(u)*cos(self.omega_big)*cos(self.i),
                         cos(u)*sin(self.i)])

    def speed_greenwich(self, minutes):
        gamma = Earth.w_sun*minutes + Earth.w_self*minutes
        return greenwich_to_equatorial(self.speed_equatorial(minutes), gamma)


class TransitSatellite(Satellite):
    def __init__(self, omega_big, tau):
        super().__init__(omega_big, tau)
        self.i = radians(90)  # degrees -> radians
        self.R = 7500  # km
        self.w = radians(3)  # degrees/min -> radians/min


class GPSSatellite(Satellite):
    def __init__(self, omega_big, tau):
        super().__init__(omega_big, tau)
        self.i = radians(60)  # degrees -> radians
        self.R = 15000  # km
        self.w = radians(2)  # degrees/min -> radians/min
