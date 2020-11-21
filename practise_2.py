from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy, find_force
from functions.math_functions import vector_module
from functions.system_evolution import verlet


particles, L = [[1, 0, 0], [1, 2, 0]], 0.8
particles2 = particles.copy()

print(particles)
for i in range(480):
    particles, particles2, vel = verlet(particles, particles2, 0.1, L, L / 2)

print(particles2)