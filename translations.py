import numpy as np
import math

def three_to_one(three_coords, gamma):
    matrix = np.array([[math.cos(gamma), -math.sin(gamma), 0],
                       [math.sin(gamma),  math.cos(gamma), 0],
                       [              0,                0, 1]])
    return matrix.dot(three_coords)


def one_to_three(one_coords, gamma):
    matrix = np.array([[ math.cos(gamma), math.sin(gamma), 0],
                       [-math.sin(gamma), math.cos(gamma), 0],
                       [               0,               0, 1]])
    return matrix.dot(one_coords)

def two_to_one(two_coords):
    R, psi, lam = two_coords
    return np.array([R*math.cos(psi)*math.cos(lam),
                     R*math.cos(psi)*math.sin(lam),
                     R*math.sin(psi)])

def one_to_two(one_coords):
    x, y, z = one_coords
    R = np.linalg.norm([x, y, z])
    psi = math.atan2(z, R)
    lam = math.atan2(y, x)
    return np.array([R, psi, lam])