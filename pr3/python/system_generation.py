import random as rng
import numpy as np


def sc_lattice_velocity(N, ro):
    a = 1 / (ro ** (1 / 3))
    L = a * N
    M = N ** 3
    particles = []
    velocities = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                particles.append([x * a, y * a, z * a])
                velocities.append(
                    [
                        rng.uniform(-0.5, 0.5),
                        rng.uniform(-0.5, 0.5),
                        rng.uniform(-0.5, 0.5),
                    ]
                )
    particles = np.asanyarray(particles)
    velocities = np.asanyarray(velocities)
    velocities = velocities - np.sum(velocities) / (M * 3)
    return particles, velocities, L
