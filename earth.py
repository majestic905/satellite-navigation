from math import radians


class Earth:
    R = 6300  # km. not bad, huh?
    w_self = radians(360/24/60)  # rad/min
    w_sun = radians(360/365/24/60)  # rad/min
