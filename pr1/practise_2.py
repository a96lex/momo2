from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy, find_force
from functions.math_functions import vector_module
from functions.system_evolution import verlet

import matplotlib.pyplot as plt

particles, L = [[0.0, 1.1, 0.0], [0.0, 0, 0.0]], 1000
particles2 = particles.copy()

test1 = []
test2 = []
test3 = []
print(particles)
for i in range(40000):
    particles, particles2, vel, pot = verlet(particles, particles2, 0.001, L, L / 2)
    test1.append(particles[0][1])
    test2.append(particles[1][1])
    test3.append(pot)

itera = []
for i in range(len(test1)):
    itera.append(i)

plt.plot(itera, test1)
plt.plot(itera, test2)
plt.show()
plt.plot(itera, test3)
plt.show()