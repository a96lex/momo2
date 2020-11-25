import numpy as np
from funcions.system_evolution import verlet
import matplotlib.pyplot as plt

part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
part2 = np.copy(part1)
L = 3
dt = 0.0001


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


def verlet(part1, part2, dt, L, cutoff):
    forces, pot = find_force_LJ0(part1, cutoff, L)
    part3 = np.zeros(np.shape(part1))
    vel = np.zeros(np.shape(part1))
    for i in range(len(part1)):
        part3[i] = 2 * part2[i] - part1[i] + forces[i] * dt ** 2
        vel[i] = (part2[i] - part1[i]) / (2 * dt)
    return part2, part3, pot


test1 = []
test2 = []
it = []
for i in range(100000):
    part1, part2, pot = verlet(part1, part2, dt, L, L / 2)
    test1.append(part1[0][1])
    test2.append(part2[0][1])
    it.append(i)


plt.plot(it, test1)
plt.plot(it, test2)
plt.show()
