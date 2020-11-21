# Modelització molecular

## Molecular Dynamics Simulations

### Practise 1

<details>
<summary>Show practise 1</summary>

##### 1. Write code (preferrably a function or subroutine) to initialize the positions of particles in a sc lattice

<details>
<summary>Show solution</summary>

In order to do so, we can directly call the previous function. In order to export the array as a .xyz readable file, I constructed a function that takes an array of particles and a title and outputs a .xyz file with its contents (always )

<details>
<summary>Show code</summary>

```python
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
```

</details>

</details>

##### 2. Prepare a system of 216 particles in a sc lattice with reduced density ρ = 0.8. Visualize and generate a snapshot of the resulting configuration (call it initconf.tga).

<details>
<summary>Show solution</summary>

The following function takes an integer M and a desnsiity of particles ro as input, and returns the particle array and the value of the simulation box L.

It generates an array containing the closest values to an input integer that satisfies the dimensions of an N x N x N sc lattice structure.

It prints a warning if the input integer does not satisfy the ideal dimensionality of the box, and computes the closest appropiate value.

<details>
<summary>Show code</summary>

Write to file function

```python
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
```

Call both functions

```python
particles, L = sc_lattice(216, 0.8)
write_file(particles,"fcc_lattice")
```

</details>

I then opened the generated file (fcc_lattice216.xyz) in jmol and got a snapshot from there

<details>
<summary>Show snapshot</summary>

Write to file function

```python
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
```

Call both functions

```python
particles, L = sc_lattice(216, 0.8)
write_file(particles,"fcc_lattice")
```

</details>

</details>
</details>
