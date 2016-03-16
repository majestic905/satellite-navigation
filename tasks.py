import numpy as np
import math
from datetime import datetime, timezone, timedelta

from translations import *
from satellites import *
from utils import *
from tests import *


transit_sats = [TransitSatellite(30, 40), TransitSatellite(60, 10), TransitSatellite(120, 100),
                TransitSatellite(210, 80), TransitSatellite(270, 50), TransitSatellite(330, 110)]

gps_sats = [GPSSatellite(0, 30), GPSSatellite(45, 100), GPSSatellite(90, 70),
            GPSSatellite(135, 0), GPSSatellite(210, 150), GPSSatellite(300, 120)]


def get_visible_satellites(where, when):
    """
    :param where: [psi, lam] radians
    :param when: datetime with timezone
    :return: [[a, b], [c, d]] transits and gpses indices
    """
    # latitude and longitude
    psi, lam = where

    utc = timezone(timedelta(hours=0))
    when_utc = when.astimezone(utc)
    print("Date and time:", when_utc.isoformat(' '))

    minutes = minutes_since_spring_equinox(when_utc)
    print("Minutes since spring equinox:", minutes)
    # next two are identical as angles, mod 2*pi
    gamma = Earth.w_sun*minutes + Earth.w_self*minutes
    gamma_test = Earth.w_sun*minutes + Earth.w_self*(when_utc.hour*60 + when_utc.minute)
    print("Gamma:", math.degrees(gamma_test))

    point_cds = one_to_three(two_to_one(np.array([Earth.R, psi, lam])), gamma)
    transit_cds = map(lambda sat: sat.three_cds(minutes), transit_sats)
    gps_cds = map(lambda sat: sat.three_cds(minutes), gps_sats)

    transit_vis = [i for i, sat_cds in enumerate(transit_cds, start=1) if visible_from_ded_position(point_cds, sat_cds)]
    print('Visible Transits:', list(transit_vis), sep='\n', end='\n\n')

    gps_vis = [i for i, sat_cds in enumerate(gps_cds, start=1) if visible_from_ded_position(point_cds, sat_cds)]
    print('Visible GPSes:', list(gps_vis), sep='\n', end='\n\n')

    return [transit_vis, gps_vis]


def rho_rho2(pos_sph2, sats_cds, distances):
    """
    :param pos_sph2: [psi, lam] radians
    :param sats_cds: [[x, y, z], [x, y, z]] rectangular
    :param distances: [x, x] km
    :param return: [psi, lam] radians
    """
    for i in range(1, 11):
        pos_sph3 = np.insert(pos_sph2, 0, [Earth.R])
        pos_rect = two_to_one(pos_sph3)
        dist_approx = np.array([distance(pos_rect, sat) for sat in sats_cds])
        print(dist_approx)
        A = np.array([distance_derivative(pos_sph3, sat, dist, True) for sat, dist in zip(sats_cds, dist_approx)])
        dd = distances - dist_approx
        dq = np.linalg.solve(A, dd)
        pos_sph2 += dq

        print("Step:", i)
        print("dq: {}, norm: {}".format(dq, np.linalg.norm(dq)))
        print("Position:", np.degrees(pos_sph2), end='\n\n')

    return pos_sph2


def rho_rho3(pos_sph3, sats_cds, distances):
    """
    :param pos_sph3: [R, psi, lam] km radians
    :param sats_cds: [[x, y, z], [x, y, z], [x, y, z]] rectangular
    :param distances: [x, x, x] km
    :param return: [R, psi, lam] km radians
    """
    pos_rect = two_to_one(pos_sph3)

    for i in range(1, 11):
        dist_approx = np.array([distance(pos_rect, sat) for sat in sats_cds])
        A = np.array([distance_derivative(pos_rect, sat, dist, False) for sat, dist in zip(sats_cds, dist_approx)])
        dd = distances - dist_approx
        dq = np.linalg.solve(A, dd)
        pos_rect += dq

        pos_sph3 = one_to_two(pos_rect)
        pos_sph3_deg = np.array([pos_sph3[0], np.degrees(pos_sph3[1]), np.degrees(pos_sph3[2])])
        print("Step:", i)
        print("dq: {}, norm: {}".format(dq, np.linalg.norm(dq)))
        print("Position:", pos_sph3_deg, end='\n\n')

    pos_sph3 = one_to_two(pos_rect)
    return pos_sph3


def get_transition_matrix(position, sats_cds, distances):
    """
    :param position: [psi, lam] radians
    :param sats_cds: [[x, y, z], [x, y, z]] rectangular
    :param distances: [x, x] km
    :return: M = [[a, b], [c, d]] such that [pos1 pos2] = M*[pos0 pos1]
    """
    position1 = rho_rho2(position, sats_cds, distances[0])
    position2 = rho_rho2(position1, sats_cds, distances[1])

    S0 = np.matrix.transpose(np.array((position, position1)))
    S1 = np.matrix.transpose(np.array((position1, position2)))

    # S1 = M*S0
    M = S1.dot(np.linalg.inv(S0))

    return M


if __name__ == '__main__':
    # test_satellites_visibility()
    # test_rho_rho2()
    # test_rho_rho3()

    initial_point = np.radians([85, 20])
    unknown_point = np.radians([87, 19])

    msk_tz = timezone(timedelta(hours=3))
    utc_tz = timezone(timedelta(hours=0))
    when_msk = datetime(2015, 2, 15, 11, 10, 0, 0, msk_tz)
    when_utc = when_msk.astimezone(utc_tz)

    _, gps_indices = get_visible_satellites(unknown_point, when_utc)

    minutes = minutes_since_spring_equinox(when_utc)
    sats_cds = np.array([gps_sats[i].one_cds(minutes) for i in gps_indices])

    unknown_point_rect = two_to_one(np.insert(unknown_point, 0, [Earth.R]))
    distances = np.array([distance(unknown_point_rect, sat) for sat in sats_cds])

    updated_point = rho_rho2(initial_point, sats_cds, distances)

    print("Point we are at:", np.degrees(unknown_point))
    print("Point we have located:", np.degrees(updated_point))
