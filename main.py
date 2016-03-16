import numpy as np
import math
from datetime import datetime, timezone, timedelta

from translations import *
from satellites import *
from utils import *


transit_sats = [TransitSatellite(30, 40), TransitSatellite(60, 10), TransitSatellite(120, 100),
                TransitSatellite(210, 80), TransitSatellite(270, 50), TransitSatellite(330, 110)]

gps_sats = [GPSSatellite(0, 30), GPSSatellite(45, 100), GPSSatellite(90, 70),
            GPSSatellite(135, 0), GPSSatellite(210, 150), GPSSatellite(300, 120)]


def task_one(where, when_msk):
    # latitude and longitude
    psi, lam = map(math.radians, where)

    utc = timezone(timedelta(hours=0))
    when_utc = when_msk.astimezone(utc)
    print("Date and time:", when_utc.isoformat(' '))

    minutes = minutes_since_spring_equinox(when_utc)
    print("Minutes since spring equinox:", minutes)
    # next two are identical as angles, mod 2*pi
    gamma = Earth.w_sun*minutes + Earth.w_self*minutes
    tchaschnikova = Earth.w_sun*minutes + Earth.w_self*(when_utc.hour*60 + when_utc.minute)
    print("Gamma:", math.degrees(tchaschnikova))

    point_cds = one_to_three(two_to_one(np.array([Earth.R, psi, lam])), gamma)
    transit_cds = map(lambda sat: sat.three_cds(minutes), transit_sats)
    gps_cds = map(lambda sat: sat.three_cds(minutes), gps_sats)

    transit_vis = [i for i, sat_cds in enumerate(transit_cds, start=1) if visible_from_ded_position(point_cds, sat_cds)]
    print('Visible Transits:', list(transit_vis), sep='\n', end='\n\n')

    gps_vis = [i for i, sat_cds in enumerate(gps_cds, start=1) if visible_from_ded_position(point_cds, sat_cds)]
    print('Visible GPSes:', list(gps_vis), sep='\n', end='\n\n')

    return [transit_vis, gps_vis]


def test_task_one():
    msk = timezone(timedelta(hours=3))

    # 1
    # Какие спутники системы ТРАНЗИТ и ЖПС будут находиться 15 февраля в 11 часов 10 минут московского времени (сдвиг 3
    # часа) в зоне видимости объекта с координатами 85 градусов северной широты и 20 градусов восточной долготы?
    when_msk = datetime(2015, 2, 15, 11, 10, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[3], [2, 3]]

    # 2
    when_msk = datetime(2015, 3, 12, 17, 10, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[3, 6], [2, 3]]

    # 3
    when_msk = datetime(2015, 4, 1, 12, 30, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[5], [4, 5]]

    # 4
    # when_msk = datetime(2015, 7, 11, 23, 00, 0, 0, msk)
    # assert task_one([85, 20], when_msk) == [[3, 4], [2, 3]]

    # 5
    when_msk = datetime(2015, 1, 22, 20, 15, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[1, 5], [2, 3, 6]]

    # 6
    when_msk = datetime(2015, 6, 24, 12, 25, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[5], [4, 5]]

    # 7
    when_msk = datetime(2015, 10, 27, 18, 5, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[1], [5, 6]]

    # 8
    when_msk = datetime(2015, 11, 2, 9, 5, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[3], [5, 6]]

    # 9
    when_msk = datetime(2015, 12, 31, 19, 5, 0, 0, msk)
    assert task_one([85, 20], when_msk) == [[3], [1, 4]]

    # 10
    # when_msk = datetime(2015, 5, 7, 11, 45, 0, 0, msk)
    # assert task_one([85, 20], when_msk) == [[2], [2, 5, 6]]

    # 11
    # when_msk = datetime(2015, 8, 2, 21, 20, 0, 0, msk)
    # assert task_one([85, 20], when_msk) == [[6], [4, 5]]

    # 12
    # when_msk = datetime(2015, 9, 1, 15, 30, 0, 0, msk)
    # assert task_one([85, 20], when_msk) == [[2], [4, 5]]


def task_two(coords, satellites, dist_real):
    coords = np.array(list(map(math.radians, coords)), dtype='float64')
    dist_real = np.array(dist_real)
    converges = False

    for i in range(1,11):
        pos_approx = np.insert(coords, 0, [Earth.R])
        dist_approx = np.array([distance(pos_approx, sat, spherical=True) for sat in satellites])
        A = np.array([distance_derivative_sph(pos_approx, sat, dist) for sat, dist in zip(satellites, dist_approx)])
        dd = dist_real - dist_approx
        dq = np.linalg.solve(A, dd)
        dq_norm = np.linalg.norm(dq)

        coords += dq
        print("Step:", i)
        print("dq: {}, norm: {}".format(dq, dq_norm))
        print("Position:", list(map(math.degrees, coords)), end='\n\n')

        if dq_norm < 1e-6:
            converges = True

    return converges

def test_task_two():
    # 1
    # Объект находится в Гвинейском заливе в окрестности точки с коорднатами 0.5 градуса северной широты и 1 градус
    # западной долготы. В зоне видимости 2 спутика с прямоугольными гринвичскими координатами x, y, z = [8676.21,
    # -2487.86, 4305.11] и x, y, z = [9146.41, 2280.46, -3338.07]. Используется дальномерный метод навигации.
    # Измеренные расстояния равны соответственно 5339.75 и 5101.44. Уточнить координаты методом Ньютона, используя
    # три итерации.
    pos_approx = [0.5, -1]
    satellites = [[8676.21, -2487.86, 4305.11],
                  [9146.41, 2280.46, -3338.07]]
    dist_real = [5339.75, 5101.44]
    assert task_two(pos_approx, satellites, dist_real) == True

    # 2
    pos_approx = [0.5, 0.5]
    satellites = [[8738.36, -2259.89, 4305.11],
                  [9083.58, 2519.10, -3338.07]]
    dist_real = [5468.63, 4986.28]
    assert task_two(pos_approx, satellites, dist_real) == True

    # 3
    pos_approx = [1, 1]
    satellites = [[8720.96, -2174.38, 4383.71]
                  [9088.91, 2606.2, -3255.68]]
    dist_real = [5196.41, 5277.57]
    assert task_two(pos_approx, satellites, dist_real) == True


def task_three(pos_approx, satellites, dist_real):
    converges = False
    dist_real = np.array(dist_real)
    pos_approx = two_to_one(np.array([pos_approx[0], math.radians(pos_approx[1]), math.radians(pos_approx[2])],
                                     dtype='float64'))

    for i in range(1, 11):
        dist_approx = np.array([distance(pos_approx, satellite, spherical=False) for satellite in satellites])
        A = np.array([distance_derivative_rect(pos_approx, sat, dist) for sat, dist in zip(satellites, dist_approx)])
        dd = dist_real - dist_approx
        dq = np.linalg.solve(A, dd)
        dq_norm = np.linalg.norm(dq)

        pos_approx += dq
        pos_approx_two = one_to_two(pos_approx)
        pos_approx_two = [pos_approx_two[0], math.degrees(pos_approx_two[1]), math.degrees(pos_approx_two[2])]
        print("Step:", i)
        print("dq: {}, norm: {}".format(dq, dq_norm))
        print("Position:", pos_approx_two, end='\n\n')

        if dq_norm < 1e-6:
            converges = True

    return converges

def test_task_three():
    # ?
    pos_approx = [7000, 1, -1]
    satellites = [[8639.76, -2477.41, 4383.71],
                  [9174.33, 2287.42, -3255.68],
                  [8974.15, -3444.85, 2756.37]]
    dist_real = [5141.61, 4691.28, 4656.23]
    task_three(pos_approx, satellites, dist_real)

if __name__ == '__main__':
    # test_task_one()
    # test_task_two()
    # test_task_three()
    pass