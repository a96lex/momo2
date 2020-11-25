from .math_functions import suma_vectors, scalar_prod_vectors
import numpy as np


def distance(p1, p2, L):
    dist = 0
    dr = []
    if type(p1) == list:
        for i in range(len(p1)):
            di = p1[i] - p2[i]
            if di > L / 2:
                dr.append(di - L)
            elif di < -L / 2:
                dr.append(di + L)
            else:
                dr.append(di)
        for i in range(len(p1)):
            dist += (dr[i]) ** 2
        dist = dist ** (1 / 2)
    else:
        di = p1 - p2
        print("d", di, L)
        if di > L / 2:
            dist = di - L
        elif di < -L / 2:
            dist = di + L
        else:
            dist = di
        print("dist", dist)
    return dist


def calculate_energy(particles, cutoff, L):
    energy = []
    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            dist = distance(particles[i], particles[j], L)
            if dist < cutoff:
                energy.append(4 * (1 / dist ** 12 - 1 / dist ** 6))

    return energy


def find_force(particles, cutoff, L):
    force_vector = []
    pot = 0
    for i in range(len(particles)):
        forcei = [0, 0, 0]
        for j in range(len(particles)):
            if i != j:
                dx = distance(particles[i][0], particles[j][0], L)
                dy = distance(particles[i][1], particles[j][1], L)
                dz = distance(particles[i][2], particles[j][2], L)
                dist = distance(particles[i], particles[j], L)
                print(dx, dy, dz)
                if dist < cutoff:
                    placeholder_dist = 48 / (dist ** 14) - 24 / (dist ** 8)
                    forceij = [
                        placeholder_dist * dx,
                        placeholder_dist * dy,
                        placeholder_dist * dz,
                    ]
                    forcei = suma_vectors(forcei, forceij)
        force_vector.append(forcei)
        pot = (
            pot
            + 4.0 * (1 / dist ** 12 - 1 / dist ** 6)
            - 4.0 * (1 / cutoff ** 12 - 1 / cutoff ** 6)
        )
    return force_vector, pot


def find_force_LJ0(r, L_box, cutoff):
    N = len(r)
    F = np.zeros((N, 3))
    pot = 0.0
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = r[i][0] - r[j][0]
                dy = r[i][1] - r[j][1]
                dz = r[i][2] - r[j][2]
                d = (
                    dx ** 2 + dy ** 2 + dz ** 2
                ) ** 0.5  # Calculating distance between particles i and j
                if d < cutoff:
                    F[i][0] = F[i][0] + (48 / d ** 14 - 24 / d ** 8) * dx
                    F[i][1] = F[i][1] + (48 / d ** 14 - 24 / d ** 8) * dy
                    F[i][2] = F[i][2] + (48 / d ** 14 - 24 / d ** 8) * dz
        pot = (
            pot
            + 2.0 * (1 / d ** 12 - 1 / d ** 6)
            - 2.0 * (1 / cutoff ** 12 - 1 / cutoff ** 6)
        )
    print("f", np.array(F))
    return np.array(F), pot