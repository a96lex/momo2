import numpy as np
from funcions.system_evolution import verlet
import matplotlib.pyplot as plt

part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
part2 = np.copy(part1)
L = 10
dt = 0.001
iteracions = 10000
m = 1


def find_force_LJ0(r, L_box, cutoff):
    N = len(r)
    F = np.zeros((N, 3))
    pot = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            dx = r[i][0] - r[j][0]
            dy = r[i][1] - r[j][1]
            dz = r[i][2] - r[j][2]
            d = (
                dx ** 2 + dy ** 2 + dz ** 2
            ) ** 0.5  # Calculating distance between particles i and j
            if d < cutoff:
                F[i][0] = F[i][0] + (48 / d ** 14 - 24 / d ** 8) * dx
                F[i][1] = F[i][1] + (48 / d ** 14 - 24 / d ** 8) * dy
                F[i][2] = F[i][2] + (48 / d ** 14 - 24 / d ** 8) * dz

                F[j][0] = F[j][0] - (48 / d ** 14 - 24 / d ** 8) * dx
                F[j][1] = F[j][1] - (48 / d ** 14 - 24 / d ** 8) * dy
                F[j][2] = F[j][2] - (48 / d ** 14 - 24 / d ** 8) * dz
        pot = (
            pot
            + 4.0 * (1 / d ** 12 - 1 / d ** 6)
            - 4.0 * (1 / cutoff ** 12 - 1 / cutoff ** 6)
        )
    return np.array(F), pot


def euler(part, vel1, dt, L, cutoff, m):
    forces, pot = find_force_LJ0(part1, cutoff, L)
    part2 = part + vel1 * dt + 0.5 * forces * dt ** 2 / m
    vel2 = vel1 + forces * dt / m
    return part2, vel2, pot


def verlet(part1, part2, dt, L, cutoff, m):
    forces, pot = find_force_LJ0(part1, cutoff, L)
    part3 = 2 * part2 - part1 + forces * dt ** 2 / m
    return part2, part3, pot


def vel_verlet(part, vel1, dt, L, cutoff, m):
    forces, pot = find_force_LJ0(part1, cutoff, L)
    part2 = part + vel1 * dt + 0.5 * forces * dt ** 2 / m
    vel2 = vel1 + forces * 0.5 * dt
    return part2, vel2, pot


test1 = []
test2 = []
test3 = []
t = np.arange(0, iteracions) * dt
vel = np.zeros(np.shape(part1))
pot1 = []
pot2 = []
pot3 = []
for i in range(iteracions):
    part1, vel, pot = euler(part1, vel, dt, L, L / 2, m)
    test1.append(abs(part1[0][1] - part1[1][1]))
    pot1.append(pot)

part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
vel = np.zeros(np.shape(part1))
for i in range(iteracions):
    part1, vel, pot = vel_verlet(part1, vel, dt, L, L / 2, m)
    test2.append(abs(part1[0][1] - part1[1][1]))
    pot2.append(pot)

part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
for i in range(iteracions):
    part1, part2, pot = vel_verlet(part1, part2, dt, L, L / 2, m)
    test3.append(abs(part1[0][1] - part1[1][1]))
    pot3.append(pot)


plt.plot(t, test1, label="euler")
plt.plot(t, test2, label="verlet vel")
plt.plot(t, test3, label="verlet")
plt.legend()
plt.show()
