from functions.system_generation import sc_lattice_velocity
from functions.system_functions import (
    vel_verlet,
    therm_Andersen,
    distance,
    time_step_vVerlet_temp,
)
from functions.write_xyz import write_file
import numpy as np


# ex 1
particles, velocities, L = sc_lattice_velocity(4, 0.8442)

# ex 2
temp = 10
nu = 1
dt = 0.1

# ex 3
particles, velocities, L = sc_lattice_velocity(3, 0.8442)

velocities = np.zeros(particles.shape)
T = 1

for i in range(20):
    write_file(particles, "test2")
    # particles, velocities, pot, kin = vel_verlet(
    #     particles, velocities, dt, L, L / 2, 1, T
    # )
    particles, velocities, pot, kin = time_step_vVerlet_temp(
        particles, 1, dt, L, T ** 0.5
    )
