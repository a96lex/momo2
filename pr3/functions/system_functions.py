import numpy as np


def therm_Andersen(velocities, nu, sigma):
    vel = np.copy(velocities)
    c = 0
    for n in range(len(vel)):
        a = np.random.ranf()
        if a < 0.1:
            c += 1
            vel[n] = np.random.normal(0, sigma, 3)
    return vel, c
