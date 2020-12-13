from functions.system_generation import sc_lattice_velocity
from functions.system_functions import therm_Andersen

particles, velocities, L = sc_lattice_velocity(2, 0.8442)

temp = 10
nu = 1

vel, c = therm_Andersen(velocities, nu, temp ** 0.5)
