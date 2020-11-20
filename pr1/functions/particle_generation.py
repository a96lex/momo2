def sc_lattice(M, ro):
    N = int(
        round(M ** (1 / 3)) // 1
    )  # Nodes of the 3d cube for n particles. As it is an integer value, the rest of the particles won't be taken care of
    res = M - N ** 3
    if res != 0:
        print(
            "The number of particles does not correspond to a N x N x N cube\nThe simulation will use the closest integer value ("
            + str(N)
            + " x "
            + str(N)
            + " x "
            + str(N)
            + ")"
        )
    a = 1 / (ro ** (1 / 3))
    L = a * N
    particles = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                particles.append([x * a, y * a, z * a])
    return particles, L


def fcc_lattice(M, ro):
    N = int(
        round((M / 4) ** (1 / 3)) // 1
    )  # Nodes of the 3d cube for n particles. As it is an integer value, the rest of the particles won't be taken care of
    res = M - N ** 3 * 4
    if res != 0:
        print(
            "The number of particles does not correspond to a 4*N^3 fcc lattice\nThe simulation will use the closest integer value (4*"
            + str(N)
            + "^3)"
        )
    a = (4 / ro) ** (1 / 3)
    L = a * N
    particles = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                particles.append([x * a, y * a, z * a])
                particles.append([x * a, y * a + a / 2, z * a + a / 2])
                particles.append([x * a + a / 2, y * a + a / 2, z * a])
                particles.append([x * a + a / 2, y * a, z * a + a / 2])
    return particles, L


def write_file(particles, filetype):
    n = len(particles)
    f = open(filetype + str(n) + ".xyz", "w")
    f.write(str(n) + "\n")
    f.write(filetype + str(n) + ".xyz\n")
    for particle in particles:
        string = ""
        for j in particle:
            string += str(j) + " "
        f.write("C " + string + "\n")
    f.close
