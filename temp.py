import numpy as np
from math import radians, degrees
from earth import *
from translations import *

if __name__ == '__main__':
    ipoint_sph = np.array([Earth.R, radians(87), radians(19)])
    rpoint_sph = np.array([Earth.R, radians(85), radians(20)])
    ipoint_rect = spherical_to_rectangular(ipoint_sph)
    rpoint_rect = spherical_to_rectangular(rpoint_sph)
    ipoint = rectangular_to_spherical(ipoint_rect)
    rpoint = rectangular_to_spherical(rpoint_rect)
    print(np.degrees(ipoint[1:3]), np.degrees(rpoint[1:3]))
