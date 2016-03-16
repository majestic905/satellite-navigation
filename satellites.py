import numpy as np
import math

from utils import *
from translations import *
from earth import *


class Satellite:
    # eps, omega = 0, 0

    def __init__(self, omega_big, tau):
        self.omega_big = math.radians(omega_big)
        self.tau = tau
        self.i, self.R, self.w = 0, 0, 0

    def one_cds(self, minutes):
        # u = self.omega + self.w*(t-self.tau)  # omega = 0
        u = self.w*(minutes - self.tau)

        in_orbit_plane = self.R*np.array([math.cos(u), math.sin(u), 0])

        around_x = np.array([[1,                0,                 0],
                             [0, math.cos(self.i), -math.sin(self.i)],
                             [0, math.sin(self.i),  math.cos(self.i)]])

        around_z = np.array([[math.cos(self.omega_big), -math.sin(self.omega_big), 0],
                             [math.sin(self.omega_big),  math.cos(self.omega_big), 0],
                             [                       0,                         0, 1]])

        return around_z.dot(around_x).dot(in_orbit_plane)

    def three_cds(self, minutes):
        sidereal_time = Earth.w_sun*minutes
        gamma = sidereal_time + Earth.w_self*minutes

        return three_to_one(self.one_cds(minutes), gamma)


class TransitSatellite(Satellite):
    def __init__(self, omega_big, tau):
        super().__init__(omega_big, tau)
        self.i = math.radians(90)  # degrees -> radians
        self.R = 7500  # km
        self.w = math.radians(3)  # degrees/min -> radians/min


class GPSSatellite(Satellite):
    def __init__(self, omega_big, tau):
        super().__init__(omega_big, tau)
        self.i = math.radians(60)  # degrees -> radians
        self.R = 15000  # km
        self.w = math.radians(2)  # degrees/min -> radians/min
