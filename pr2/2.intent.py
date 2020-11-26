import numpy as np
from funcions.write_xyz import write_file
from funcions.system_functions import find_force_LJ0, vel_verlet, verlet, euler
import matplotlib.pyplot as plt

part1 = np.array([[1, 1.1, 0], [0, 0, 0]])
part2 = np.copy(part1)
L = 10
dt = 0.0003
iteracions = 30000
m = 1


pose = []
posv = []
posvv = []
pote = []
potv = []
potvv = []
cine = []
cinv = []
cinvv = []
enee = []
enev = []
enevv = []


t = np.arange(0, iteracions) * dt
vel = np.zeros(np.shape(part1))
for i in range(iteracions):
    # write_file(part1, "teest")
    part1, vel, pot, kin = euler(part1, vel, dt, L, L / 2, m)
    pose.append(abs(part1[0][1] - part1[1][1]))
    enee.append(pot + kin)
    pote.append(pot)
    cine.append(kin)


part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
vel = np.zeros(np.shape(part1))
for i in range(iteracions):
    part1, vel, pot, kin = vel_verlet(part1, vel, dt, L, L / 2, m)
    posvv.append(abs(part1[0][1] - part1[1][1]))
    enevv.append(pot + kin)
    potvv.append(pot)
    cinvv.append(kin)


part1 = np.array([[0, 1.1, 0], [0, 0, 0]])
part2 = np.array([[0, 1.1, 0], [0, 0, 0]])
for i in range(iteracions):
    part1, part2, pot, kin = verlet(part1, part2, dt, L, L / 2, m)
    posv.append(abs(part1[0][1] - part1[1][1]))
    enev.append(pot + kin)
    potv.append(pot)
    cinv.append(kin)


plt.plot(t, enee, label="euler")
plt.plot(t, enev, label="verlet")
plt.plot(t, enevv, label="verlet vel")
plt.legend()
plt.show()
