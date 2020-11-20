from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy, find_force
from functions.math_functions import vector_module


particles, L = sc_lattice(216, 0.8)
energy_vector = calculate_energy(particles, L / 2, L)
energy_matrix = find_force(particles, L / 2, L)
print("sc lattice for " + str(len(particles)) + " particles")
print("Lennard-Jones potential energy (with perdiodic boundary conditions)")
print("cutoff: ", L / 2, " energy: ", vector_module(energy_vector))
write_file(particles, "test")
