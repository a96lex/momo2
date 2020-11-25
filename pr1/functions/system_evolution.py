from .system_functions import find_force, find_force_LJ0
from .math_functions import suma_vectors, scalar_prod_vectors, suma_vectors_prod


def verlet(part1, part2, dt, L, cutoff):
    forces, pot = find_force(part1, cutoff, L)
    print("f", forces)
    part3 = []
    vel = []
    for i in range(len(part1)):
        vel.append(
            scalar_prod_vectors(
                suma_vectors(part1[i], scalar_prod_vectors(part2[i], -1)), 1 / (2 * dt)
            )
        )

        # (part2[i]*2+part1[i]*-1)*1+forces[i]* dt ** 2 ,
        part3.append(
            suma_vectors_prod(
                suma_vectors_prod(part2[i], 2, part1[i], -1),
                1,
                forces[i],
                dt ** 2,
            )
        )
    return part2, part3, vel, pot
