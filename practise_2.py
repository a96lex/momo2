from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy, find_force
from functions.math_functions import vector_module
from functions.system_evolution import verlet

import matplotlib.pyplot as plt

particles, L = [[0.0, 1.283879891646726, 0.0], [0.0, -0.28387989164675953, 0.0]], 3
particles2 = particles.copy()

test1 = []
test2 = []
test3 = []
print(particles)
for i in range(10000):
    particles, particles2, vel, pot = verlet(particles, particles2, 0.001, L, L / 2)
    print(i, particles, "\n")
    test1.append(particles[0][1])
    test2.append(particles[1][1])

itera = []
for i in range(len(test1)):
    itera.append(i)

plt.plot(itera, test1)
plt.plot(itera, test2)
plt.show()
# plt.plot(itera, test3)
# plt.show()