from functions.particle_generation import sc_lattice, fcc_lattice, write_file

particles, L = sc_lattice(216, 0.8)
write_file(particles, "xyz/sc_lattice_")
particles, L = fcc_lattice(256, 0.8)
write_file(particles, "xyz/fcc_lattice_")
