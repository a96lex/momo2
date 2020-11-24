from .math_functions import suma_vectors


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
    else:
        di = p1 - p2
        if di > L / 2:
            return di - L
        elif di < -L / 2:
            return di + L
        else:
            return di
    dist = dist ** (1 / 2)
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
    for i in range(len(particles)):
        forcei = [0, 0, 0]
        for j in range(len(particles)):
            if i != j:
                dx = distance(particles[i][0], particles[j][0], L)
                dy = distance(particles[i][1], particles[j][1], L)
                dz = distance(particles[i][2], particles[j][2], L)
                dist = distance(particles[i], particles[j], L)
                if dist < cutoff:
                    placeholder_dist = 48 / dist ** 14 - 24 / dist ** 8
                    forceij = [
                        placeholder_dist * dx,
                        placeholder_dist * dy,
                        placeholder_dist * dz,
                    ]
                    forcei = suma_vectors(forcei, forceij)
            force_vector.append(forcei)
    return force_vector