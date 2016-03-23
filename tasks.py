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
    :param display: bool - defines whether to display information or not
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
    :param display: bool - defines whether to display information or not
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


def rho_rho3(position, sats_cds, distances, display=False):
    """
    :param position: [R, psi, lam] km radians - approximate position, R isn't necessary Earth's radius
    :param sats_cds: [[x, y, z], [x, y, z], [x, y, z]] rectangular - visible satellites coordinates
    :param distances: [x, x, x] km - measured distances to satellites from position
    :param display: bool - defines whether to display information or not
    :return: [R, psi, lam] km radians - more accurate position
    """
    position = spherical_to_rectangular(position)

    for i in range(1, 11):
        dist_approx = np.array([distance(position, sat) for sat in sats_cds])
        A = np.array([distance_derivative(position, sat, dist, False) for sat, dist in zip(sats_cds, dist_approx)])
        dd = distances - dist_approx
        dq = np.linalg.solve(A, dd)
        position += dq

        if display:
            pos_sph3_rad = rectangular_to_spherical(position)
            pos_sph3_deg = np.array([pos_sph3_rad[0], degrees(pos_sph3_rad[1]), degrees(pos_sph3_rad[2])])
            print("Step:", i)
            print("dq: {}, norm: {}".format(dq, np.linalg.norm(dq)))
            print("Position:", pos_sph3_deg, end='\n\n')

    return rectangular_to_spherical(position)


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


def doppler(position, speed, sats_positions, sats_speeds, doppler_value_measured, display=False):
    """
    :param position: [R, psi, lam] km radians
    :param speed: [xdot, ydot, zdot] ?
    :param sats_positions: [[x, y, z], [x, y, z], [x, y, z]] - in greenwich coordinates
    :param sats_speeds: [[xdot, ydot, zdot], ...] - in greenwich coordinates
    :param doppler_value_measured: [rhodot1, rhodot2, rhodot3] - for each satellite
    :param display: bool - defines whether to display information or not
    :return: [R, psi, lam] km radians - more accurate position
    """
    position = spherical_to_rectangular(position)

    for i in range(1, 11):
        sats_dists = [distance(position, sat_pos) for sat_pos in sats_positions]
        doppler_value_computed = np.array([rho_dot(position, sat_pos, speed, sat_speed, sat_dist)
                                for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)])
        drho_dot = doppler_value_measured - doppler_value_computed
        B = np.array([rho_dot_derivative(position, sat_pos, speed, sat_speed, sat_dist)
                        for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)])
        dq = np.linalg.solve(B, drho_dot)
        dq_norm = np.linalg.norm(dq)
        position += dq

        if display:
            print("Step:", i)
            print("dq: {}, norm: {}".format(dq, dq_norm))
            print("Position:", np.degrees(rectangular_to_spherical(position)[1:3]), end='\n\n')

        if dq_norm < 1e-6:
            break

    return rectangular_to_spherical(position)


if __name__ == '__main__':
    # test_satellites_visibility()
    # test_rho_rho2()
    # test_rho_rho3()
    # test_transition_matrix()
    test_doppler()
