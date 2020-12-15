from functions.system_generation import sc_lattice_velocity
from functions.system_functions import (
    therm_Andersen,
    distance,
    time_step_vVerlet_temp,
)
from functions.write_xyz import write_file
import numpy as np
import matplotlib.pyplot as plt


# ex 1
particles, velocities, L = sc_lattice_velocity(8, 0.8442)

# ex 2
temp = 10
nu = 1
dt = 0.0001

# ex 3
# Initial conditios
lattice = 4
particles, velocities, L = sc_lattice_velocity(lattice, 0.8442)
T = 1000
n = 100
velocities = 0
temps = [0]
temp = [T]

# Setting the storage of the results
filename = "100_2.dat"
write_file(particles, filename, "w")

# Main loop
for i in range(n):
    print("\n\n\n\n\n\n\n\n\nSimulation at", 100 * i / n, "%")
    temps.append(i * dt)
    particles, velocities, kin = time_step_vVerlet_temp(
        particles, velocities, dt, L, T ** 0.5
    )
    write_file(particles, filename, "a")
    temp.append((kin * 2 / (3 * lattice ** 3)))
print("\n\n\n\n\n\n\n\n\nSimulation complete")

# Plot the temperature over time
plt.plot(temps, temp)
plt.title("Evolution of temperature over time")
plt.xlabel("Time [s]")
plt.ylabel("Reduced temperature")
plt.show()