from functions.system_generation import sc_lattice_velocity
from functions.system_functions import (
    therm_Andersen,
    distance,
    time_step_vVerlet_temp,
)
from functions.write_xyz import write_file
import numpy as np
import matplotlib.pyplot as plt


# Initial conditios
lattice = 3
dt = 0.0001
T = 1000
n = 100

# System generation
particles, velocities, L = sc_lattice_velocity(lattice, 0.8442)
velocities = 0
temps = [0]
temp = [T]
# Setting the storage of the results
filename = "100_2.dat"
write_file(particles, filename, "w")
# Main loop
for i in range(n):
    # print("\n\n\n\n\n\n\n\n\nSimulation at", 100 * i / n, "%")
    temps.append(i * dt)
    particles, velocities, T = time_step_vVerlet_temp(
        particles, velocities, dt, L, temp[0], lattice ** 3
    )
    write_file(particles, filename, "a")
    temp.append(T)
# print("\n\n\n\n\n\n\n\n\nSimulation complete")

# Plot the temperature over time
plt.plot(temps, temp)
plt.title("Evolution of temperature over time")
plt.xlabel("Time [s]")
plt.ylabel("Reduced temperature")
plt.show()