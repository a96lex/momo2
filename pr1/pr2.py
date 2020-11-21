from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy, find_force
from functions.math_functions import vector_module
from functions.system_evolution import verlet


particles, L = sc_lattice(216, 0.8)
energy_matrix = find_force(particles, L / 2, L)
print(energy_matrix)

particles, particles2, vel = verlet(particles, particles, 0.1, L, L / 2)
energy = calculate_energy(particles, L / 2, L)
