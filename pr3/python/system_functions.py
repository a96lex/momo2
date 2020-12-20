import numpy as np


def therm_Andersen(vel, nu, sigma):
    for n in range(len(vel)):
        if np.random.ranf() < nu:
            vel[n] = np.random.normal(0, sigma, 3)
    return vel


def distance(p1, p2, L):
    x = p1 - p2
    d = 0
    for i in range(len(x)):
        if x[i] > L / 2.0:
            x[i] -= L
        if x[i] < -L / 2.0:
            x[i] += L
        d += x[i] ** 2
    d = d ** 0.5
    return x, d


def find_force_LJ0(r, L, cutoff):
    N = len(r)
    F = np.zeros((N, 3))
    pot = 0.0
    c = 0
    for i in range(N):
        for j in range(i + 1, N):
            dr, d = distance(r[i], r[j], L)
            if d < cutoff:
                for k in range(3):
                    F[i][k] += (48 / d ** 14 - 24 / d ** 8) * dr[k]
                    F[j][k] -= (48 / d ** 14 - 24 / d ** 8) * dr[k]
            pot = (
                pot
                + 4.0 * (1 / d ** 12 - 1 / d ** 6)
                - 4.0 * (1 / cutoff ** 12 - 1 / cutoff ** 6)
            )
    return np.array(F), pot


def kinetic_energy(vel):
    kin = 0.0
    for i in range(len(vel)):
        kin += 0.5 * (vel[i][0] ** 2 + vel[i][1] ** 2 + vel[i][2] ** 2)
    return kin


def time_step_vVerlet_temp(r, vel, dt, L, T, N):
    F, pot = find_force_LJ0(r, L, L / 2)
    r = r + vel * dt + 0.5 * F * dt ** 2.0
    vel = vel + F * 0.5 * dt
    F, pot = find_force_LJ0(r, L, L / 2)
    vel = vel + F * 0.5 * dt
    vel = therm_Andersen(vel, 0.1, T ** 0.5)
    kin = kinetic_energy(vel)
    T = kin * 2 / (3 * N - 3)
    return r, vel, T, pot, kin
