import numpy as np
from math import sin, cos, radians, degrees
from datetime import datetime, timezone, timedelta

from translations import *
from satellites import *
from utils import *
from tests import *


transit_sats = [TransitSatellite(30, 40), TransitSatellite(60, 10), TransitSatellite(120, 100),
                TransitSatellite(210, 80), TransitSatellite(270, 50), TransitSatellite(330, 110)]

gps_sats = [GPSSatellite(0, 30), GPSSatellite(45, 100), GPSSatellite(90, 70),
            GPSSatellite(135, 0), GPSSatellite(210, 150), GPSSatellite(300, 120)]


def get_visible_satellites(where, when, display=False):
    """
    :param where: [psi, lam] radians - position spherical coordinates
    :param when: datetime.datetime with timezone - date and time
    :return: [[a, b], [c, d]] - Transit and GPS satellites indices
    """
    # TODO: what are thees? latitude and longitude?
    psi, lam = where

    utc = timezone(timedelta(hours=0))
    when_utc = when.astimezone(utc)
    minutes = minutes_since_spring_equinox(when_utc)
    if display:
        print("Date and time:", when_utc.isoformat(' '))
        print("Minutes since spring equinox:", minutes)

    # next two are identical as angles, mod 2*pi
    gamma = Earth.w_sun*minutes + Earth.w_self*minutes
    gamma_test = Earth.w_sun*minutes + Earth.w_self*(when_utc.hour*60 + when_utc.minute)
    if display:
        print("Gamma:", degrees(gamma_test))

    point_cds = spherical_to_rectangular(np.array([Earth.R, psi, lam]))
    transit_cds = [sat.coordinates_greenwich(minutes) for sat in transit_sats]
    gps_cds = [sat.coordinates_greenwich(minutes) for sat in gps_sats]

    transit_vis = [i+1 for i, sat_cds in enumerate(transit_cds) if sat_really_visible(point_cds, sat_cds)]
    gps_vis = [i+1 for i, sat_cds in enumerate(gps_cds) if sat_really_visible(point_cds, sat_cds)]
    if display:
        print('Visible Transits:', list(transit_vis), sep='\n', end='\n\n')
        print('Visible GPSes:', list(gps_vis),  sep='\n', end='\n\n')

    return [transit_vis, gps_vis]


def rho_rho2(pos_sph2, sats_cds, distances, display=False):
    """
    :param pos_sph2: [psi, lam] radians - approximate position
    :param sats_cds: [[x, y, z], [x, y, z]] rectangular - visible satellites coordinates
    :param distances: [x, x] km - measured distances to satellites from position
    :return: [psi, lam] radians - more accurate position
    """
    for i in range(1, 11):
        pos_sph3 = np.insert(pos_sph2, 0, [Earth.R])
        pos_rect = spherical_to_rectangular(pos_sph3)
        dist_approx = np.array([distance(pos_rect, sat) for sat in sats_cds])
        A = np.array([distance_derivative(pos_sph3, sat, dist, True) for sat, dist in zip(sats_cds, dist_approx)])
        dd = distances - dist_approx
        dq = np.linalg.solve(A, dd)
        pos_sph2 += dq

        if display:
            print("Step:", i)
            print("dq: {}, norm: {}".format(dq, np.linalg.norm(dq)))
            print("Position:", np.degrees(pos_sph2), end='\n\n')

    return pos_sph2


def rho_rho3(pos_sph3, sats_cds, distances, display=False):
    """
    :param pos_sph3: [R, psi, lam] km radians - approximate position, R isn't necessary Earth's radius
    :param sats_cds: [[x, y, z], [x, y, z], [x, y, z]] rectangular - visible satellites coordinates
    :param distances: [x, x, x] km - measured distances to satellites from position
    :return: [R, psi, lam] km radians - more accurate position
    """
    pos_rect = spherical_to_rectangular(pos_sph3)

    for i in range(1, 11):
        dist_approx = np.array([distance(pos_rect, sat) for sat in sats_cds])
        A = np.array([distance_derivative(pos_rect, sat, dist, False) for sat, dist in zip(sats_cds, dist_approx)])
        dd = distances - dist_approx
        dq = np.linalg.solve(A, dd)
        pos_rect += dq

        if display:
            pos_sph3 = rectangular_to_spherical(pos_rect)
            pos_sph3_deg = np.array([pos_sph3[0], np.degrees(pos_sph3[1]), np.degrees(pos_sph3[2])])
            print("Step:", i)
            print("dq: {}, norm: {}".format(dq, np.linalg.norm(dq)))
            print("Position:", pos_sph3_deg, end='\n\n')

    pos_sph3 = rectangular_to_spherical(pos_rect)
    return pos_sph3


def get_transition_matrix(position, sats_cds, distances):
    """
    :param position: [psi, lam] radians
    :param sats_cds: [[x, y, z], [x, y, z]] rectangular
    :param distances: [x, x] km
    :return: M = [[a, b], [c, d]] such that [pos1 pos2] = M*[pos0 pos1]
    """
    position1 = rho_rho2(np.copy(position), sats_cds, distances[0], display=False)
    position2 = rho_rho2(np.copy(position1), sats_cds, distances[1], display=False)

    S0 = np.matrix.transpose(np.array((position, position1)))
    S1 = np.matrix.transpose(np.array((position1, position2)))

    # S1 = M*S0
    M = S1.dot(np.linalg.inv(S0))

    return M


def doppler():
    ipoint_sph = np.array([Earth.R, radians(87), radians(19)])
    rpoint_sph = np.array([Earth.R, radians(85), radians(20)])
    ipoint_rect = spherical_to_rectangular(ipoint_sph)
    rpoint_rect = spherical_to_rectangular(rpoint_sph)
    print("Real point:", np.degrees(rpoint_sph[1:3]), rpoint_rect)
    print("Initial point:", np.degrees(ipoint_sph[1:3]), ipoint_rect)

    msk_tz = timezone(timedelta(hours=3))
    utc_tz = timezone(timedelta(hours=0))
    when_msk = datetime(2015, 7, 11, 23, 00, 0, 0, msk_tz)
    when_utc = when_msk.astimezone(utc_tz)
    minutes = minutes_since_spring_equinox(when_utc)
    print("Date and time:", when_utc.isoformat(' '))
    print("Minutes passed since vernal equinox:", minutes)

    transit_indices, gps_indices = get_visible_satellites(rpoint_sph[1:3], when_utc)
    print("Visible satellites: Transit", transit_indices, "GPS", gps_indices)

    visible_transits = [transit_sats[i-1] for i in transit_indices]
    visible_gpses = [gps_sats[i-1] for i in gps_indices]
    selected_sats = (visible_gpses + visible_transits)[0:3]
    print("Selected satellites (GPS + Transit):", (gps_indices + transit_indices)[0:3], end='\n\n')

    sats_positions = [sat.coordinates_greenwich(minutes) for sat in selected_sats]
    sats_speeds = [sat.speed_greenwich(minutes) for sat in selected_sats]
    sats_dists = [distance(rpoint_rect, sat_pos) for sat_pos in sats_positions]

    rho_dot_real = [rho_dot(rpoint_rect, sat_pos, np.array([0, 0, 0]), sat_speed, sat_dist)
                    for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)]

    for i in range(1, 11):
        sats_dists = np.array([distance(ipoint_rect, sat_position) for sat_position in sats_positions])
        rho_dot_approx = np.array([rho_dot(ipoint_rect, sat_pos, np.array([0, 0, 0]), sat_speed, sat_dist)
                                for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)])
        drho_dot = rho_dot_real - rho_dot_approx
        B = np.array([rho_dot_derivative(ipoint_rect, sat_pos, np.array([0, 0, 0]), sat_speed, sat_dist)
                        for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)])
        dq = np.linalg.solve(B, drho_dot)
        dq_norm = np.linalg.norm(dq)
        ipoint_rect += dq

        ipoint_sph2 = rectangular_to_spherical(ipoint_rect)[1:3]
        print("Step:", i)
        print("dq: {}, norm: {}".format(dq, dq_norm))
        print("Position:", np.degrees(ipoint_sph2), end='\n\n')

        if dq_norm < 1e-6:
            break


if __name__ == '__main__':
    # test_satellites_visibility()
    # test_rho_rho2()
    # test_rho_rho3()
    # test_transition_matrix()
    # TODO: test_doppler()
    doppler()

    # initial_point = np.radians([85, 20])
    # unknown_point = np.radians([87, 19])
    #
    # msk_tz = timezone(timedelta(hours=3))
    # utc_tz = timezone(timedelta(hours=0))
    # when_msk = datetime(2015, 2, 15, 11, 10, 0, 0, msk_tz)
    # when_utc = when_msk.astimezone(utc_tz)
    #
    # _, gps_indices = get_visible_satellites(unknown_point, when_utc)
    #
    # minutes = minutes_since_spring_equinox(when_utc)
    # sats_cds = np.array([gps_sats[i].one_cds(minutes) for i in gps_indices])
    #
    # unknown_point_rect = two_to_one(np.insert(unknown_point, 0, [Earth.R]))
    # distances = np.array([distance(unknown_point_rect, sat) for sat in sats_cds])
    #
    # updated_point = rho_rho2(initial_point, sats_cds, distances)
    #
    # print("Point we are at:", np.degrees(unknown_point))
    # print("Point we have located:", np.degrees(updated_point))
