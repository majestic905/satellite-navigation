from tasks import *


def test_satellites_visibility():
    msk = timezone(timedelta(hours=3))

    # 1
    # Какие спутники системы ТРАНЗИТ и ЖПС будут находиться 15 февраля в 11 часов 10 минут московского времени (сдвиг 3
    # часа) в зоне видимости объекта с координатами 85 градусов северной широты и 20 градусов восточной долготы?
    when = datetime(2015, 2, 15, 11, 10, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[3], [2, 3]]

    # 2
    when = datetime(2015, 3, 12, 17, 10, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[3, 6], [2, 3]]

    # 3
    when = datetime(2015, 4, 1, 12, 30, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[5], [4, 5]]

    # 4
    when = datetime(2015, 7, 11, 23, 00, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[3, 4], [2, 3]]

    # 5
    when = datetime(2015, 1, 22, 20, 15, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[1, 5], [2, 3, 6]]

    # 6
    when = datetime(2015, 6, 24, 12, 25, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[5], [4, 5]]

    # 7
    when = datetime(2015, 10, 27, 18, 5, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[1], [5, 6]]

    # 8
    when = datetime(2015, 11, 2, 9, 5, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[3], [5, 6]]

    # 9
    when = datetime(2015, 12, 31, 19, 5, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[3], [1, 4]]

    # 10
    when = datetime(2015, 5, 7, 11, 45, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[2], [2, 5, 6]]

    # 11
    when = datetime(2015, 8, 2, 21, 20, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[6], [4, 5]]

    # 12
    when = datetime(2015, 9, 1, 15, 30, 0, 0, msk)
    where = np.radians([85, 20])
    assert get_visible_satellites(where, when, display=True) == [[2], [4, 5]]


def test_rho_rho2():
    # 1
    # Объект находится в Гвинейском заливе в окрестности точки с коорднатами 0.5 градуса северной широты и 1 градус
    # западной долготы. В зоне видимости 2 спутика с прямоугольными гринвичскими координатами x, y, z = [8676.21,
    # -2487.86, 4305.11] и x, y, z = [9146.41, 2280.46, -3338.07]. Используется дальномерный метод навигации.
    # Измеренные расстояния равны соответственно 5339.75 и 5101.44. Уточнить координаты методом Ньютона, используя
    # три итерации.
    position = np.radians([0.5, -1])
    satellites = [[8676.21, -2487.86, 4305.11],
                  [9146.41, 2280.46, -3338.07]]
    distances = np.array([5339.75, 5101.44])
    # assert task_two(position, satellites, distances) == True
    rho_rho2(position, satellites, distances, display=True)

    # 2
    position = np.radians([0.5, 0.5])
    satellites = [[8738.36, -2259.89, 4305.11],
                  [9083.58, 2519.10, -3338.07]]
    distances = np.array([5468.63, 4986.28])
    # assert task_two(position, satellites, distances) == True
    rho_rho2(position, satellites, distances, display=True)

    # 3
    position = np.radians([1, 1])
    satellites = [[8720.96, -2174.38, 4383.71],
                  [9088.91, 2606.2, -3255.68]]
    distances = np.array([5196.41, 5277.57])
    # assert task_two(position, satellites, distances) == True
    rho_rho2(position, satellites, distances, display=True)


def test_rho_rho3():
    # ?
    position = np.array([7000, radians(1), radians(-1)])
    satellites = np.array([[8639.76, -2477.41, 4383.71],
                           [9174.33, 2287.42, -3255.68],
                           [8974.15, -3444.85, 2756.37]])
    distances = np.array([5141.61, 4691.28, 4656.23])
    rho_rho3(position, satellites, distances, display=True)


def test_transition_matrix():
    # 1
    coords = np.radians([-1, 1])
    satellites = np.array([[8864.09, -2210.07, 4067.37],
                           [8974.15, 2573.3, -3583.68]])
    dist_real = np.array([[5466.72, 4992.39],
                          [5510.15, 4955.85]])
    matrix1 = get_transition_matrix(coords, satellites, dist_real)
    print(matrix1, end='\n\n')

    # 3
    coords = np.radians([0.5, 1])
    satellites = np.array([[8757.75, -2183.55, 4305.11],
                           [9061.25, 2598.27, -3338.07]])
    dist_real = np.array([[5506.08, 4690],
                          [5594.37, 4889.03]])
    matrix3 = get_transition_matrix(coords, satellites, dist_real)
    print(matrix3, end='\n\n')

    # 8
    coords = np.radians([0.5, 0.5])
    satellites = np.array([[8738.36, -2259.89, 4305.11],
                           [9083.58, 2519.1, -3338.07]])
    dist_real = np.array([[5403.34, 5045.4],
                          [5378.9, 5065.14]])
    matrix8 = get_transition_matrix(coords, satellites, dist_real)
    print(matrix8, end='\n\n')

    # 9
    coords = np.radians([1, 0.5])
    satellites = np.array([[8701.65, -2250.4, 4383.71],
                           [9111.3, 2526.79, -3255.68]])
    dist_real = np.array([[5389.83, 5059.78],
                          [5357.92, 5088.64]])
    matrix9 = get_transition_matrix(coords, satellites, dist_real)

    print(matrix9, sep='\n\n')


def test_doppler():
    ipoint_sph = np.array([Earth.R, radians(87), radians(19)])
    rpoint_sph = np.array([Earth.R, radians(85), radians(20)])
    print("Real point:", np.degrees(rpoint_sph[1:3]))
    print("Initial point:", np.degrees(ipoint_sph[1:3]))

    msk_tz = timezone(timedelta(hours=3))
    utc_tz = timezone(timedelta(hours=0))
    when_msk = datetime(2015, 2, 15, 11, 10, 0, 0, msk_tz)
    when_utc = when_msk.astimezone(utc_tz)
    minutes = minutes_since_spring_equinox(when_utc)
    print("Date and time:", when_utc.isoformat(' '))
    print("Minutes passed since vernal equinox:", minutes)

    transit_indices, gps_indices = get_visible_satellites(rpoint_sph[1:3], when_utc)
    visible_transits = [transit_sats[i-1] for i in transit_indices]
    visible_gpses = [gps_sats[i-1] for i in gps_indices]
    selected_sats = (visible_gpses + visible_transits)[0:3]
    print("Visible satellites: Transit", transit_indices, "GPS", gps_indices)
    print("Selected satellites (GPS + Transit):", (gps_indices + transit_indices)[0:3], end='\n\n')

    rpoint_rect = spherical_to_rectangular(rpoint_sph)
    sats_positions = [sat.coordinates_greenwich(minutes) for sat in selected_sats]
    sats_speeds = [sat.speed_greenwich(minutes) for sat in selected_sats]
    sats_dists = [distance(rpoint_rect, sat_pos) for sat_pos in sats_positions]

    rho_dot_real = [rho_dot(rpoint_rect, sat_pos, np.array([0, 0, 0]), sat_speed, sat_dist)
                    for (sat_pos, sat_speed, sat_dist) in zip(sats_positions, sats_speeds, sats_dists)]

    rpoint_sph_doppler = doppler(ipoint_sph, np.array([0, 0, 0]), sats_positions, sats_speeds, rho_dot_real, display=True)
    print("Real point:", np.degrees(rpoint_sph[1:3]))
    print("Computed point:", np.degrees(rpoint_sph_doppler[1:3]))
