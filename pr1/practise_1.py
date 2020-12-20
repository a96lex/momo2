from functions.particle_generation import sc_lattice, fcc_lattice, write_file
from functions.system_functions import calculate_energy
from functions.math_functions import vector_module


particles, L = sc_lattice(64, 0.8)
print("sc lattice for " + str(len(particles)) + " particles")
print("Lennard-Jones potential energy (no perdiodic boundary conditions)")
for i in [1.5, 2, 2.5, 3]:
    print("cutoff: ", i, " energy: ", vector_module(calculate_energy(particles, i, 0)))
print("Lennard-Jones potential energy (perdiodic boundary conditions)")
for i in [1.5, 2, 2.5, 3]:
    print("cutoff: ", i, " energy: ", vector_module(calculate_energy(particles, i, L)))


particles, L = fcc_lattice(64, 0.8)
print("\nfcc lattice for " + str(len(particles)) + " particles")
print("Lennard-Jones potential energy (no perdiodic boundary conditions)")
for i in [1.5, 2, 2.5, 3]:
    print("cutoff: ", i, " energy: ", vector_module(calculate_energy(particles, i, 0)))
print("Lennard-Jones potential energy (perdiodic boundary conditions)")
for i in [1.5, 2, 2.5, 3]:
    print("cutoff: ", i, " energy: ", vector_module(calculate_energy(particles, i, L)))
write_file(particles, "fcc2")
