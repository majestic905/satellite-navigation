import numpy as np
import math

from translations import *
from satellites import *

import datetime

def get_time(when):
    month, day, hour, minute = when.month, when.day, when.hour, when.minute
    b = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]
    return ((sum(b[0:month])+day-79) % 365)*1440 + ((hour-3) % 24)*60 + minute

if __name__ == '__main__':
    start = datetime.datetime(2015, 8, 31, 15, 30, 0)
    dt = datetime.timedelta(1)

    for i in range(100):
        today = start + i*dt
        time = get_time(today)
        sidereal_time = Earth.w_sun*time
        gamma = sidereal_time + Earth.w_self*time
        print(today.day, today.month)
        print(math.degrees(gamma % 2*math.pi))
