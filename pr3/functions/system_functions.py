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
    for i in range(N):
        for j in range(i + 1, N):
            dr, d = distance(r[i], r[j], L)
            if d < cutoff:
                F[i][0] = F[i][0] + (48 / d ** 14 - 24 / d ** 8) * dr[0]
                F[i][1] = F[i][1] + (48 / d ** 14 - 24 / d ** 8) * dr[1]
                F[i][2] = F[i][2] + (48 / d ** 14 - 24 / d ** 8) * dr[2]

                F[j][0] = F[j][0] - (48 / d ** 14 - 24 / d ** 8) * dr[0]
                F[j][1] = F[j][1] - (48 / d ** 14 - 24 / d ** 8) * dr[1]
                F[j][2] = F[j][2] - (48 / d ** 14 - 24 / d ** 8) * dr[2]
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


def euler(part, vel1, dt, L, cutoff, m):
    forces, pot = find_force_LJ0(part, cutoff, L)
    part2 = part + vel1 * dt + 0.5 * forces * dt ** 2 / m
    vel2 = vel1 + forces * dt / m
    kin = kinetic_energy(vel2)
    return part2, vel2, pot, kin


def vel_verlet(part1, vel1, dt, L, cutoff, m, T):
    forces, pot = find_force_LJ0(part1, L, cutoff)
    part2 = part1 + vel1 * dt + 0.5 * forces * dt ** 2 / m
    forces2, pot = find_force_LJ0(part2, L, cutoff)
    vel2 = vel1 + (forces + forces2) * 0.5 * dt
    vel2, c = therm_Andersen(vel2, 0.1, T)
    kin = kinetic_energy(vel2)
    return part2, vel2, pot, kin


def time_step_vVerlet_temp(r, vel, dt, L, sigma):
    F, pot = find_force_LJ0(r, L, L / 2)
    r = r + vel * dt + 0.5 * F * dt ** 2
    vel = vel + F * 0.5 * dt
    F, pot = find_force_LJ0(r, L, L / 2)
    vel = vel + F * 0.5 * dt
    vel = therm_Andersen(vel, 0.1, sigma)
    return r, vel, F, pot


def verlet(part1, part2, dt, L, cutoff, m):
    forces, pot = find_force_LJ0(part2, cutoff, L)
    part3 = 2.0 * part2 - part1 + forces * dt ** 2.0 / m
    vel = (part3 - part1) / (2.0 * dt)
    kin = kinetic_energy(vel)
    return part2, part3, pot, kin