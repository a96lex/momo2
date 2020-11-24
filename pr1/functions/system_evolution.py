from .system_functions import find_force
from .math_functions import suma_vectors, scalar_prod_vectors


def verlet(part1, part2, dt, L, cutoff):
    forces = find_force(part1, L / 2, L)
    part3 = []
    vel = []
    for i in range(len(part1)):
        part3.append(
            suma_vectors(
                suma_vectors(
                    scalar_prod_vectors(part1[i], 2), scalar_prod_vectors(part2[i], -1)
                ),
                scalar_prod_vectors(forces[i], dt ** 2),
            )
        )
        vel.append(
            scalar_prod_vectors(
                suma_vectors(part1[i], scalar_prod_vectors(part2[i], -1)), 1 / (2 * dt)
            )
        )
    return part2, part3, vel
