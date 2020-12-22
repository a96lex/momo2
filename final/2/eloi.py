import numpy as np
import matplotlib.pyplot as plt

# Control variables
input_file = "Thermodynamics.dat"
N_thermo = 300  # Nº of data before stabilization (ignored)

# Initialization of some variables
time = []
temp = []
kin = []
pot = []
tot = []
N = 0

# Data reading
with open(input_file, "r") as file:
    lines = file.readlines()

# Data parsing and storing
for line in lines:
    if line.lstrip().startswith("#"):
        continue
    elif line.lstrip() == "":
        continue
    else:
        N += 1
        t, tem, k, v, e = line.split()
        time.append(float(t))
        temp.append(float(tem))
        kin.append(float(k))
        pot.append(float(v))
        tot.append(float(e))

data = np.array(
    [
        time,
        temp,
        kin,
        pot,
        tot,
    ]
).transpose()  # Guardem les dades en les columnes data = [ time | temp | kin | pot | tot ]

# Data processing
N_true = N - N_thermo
m = 1
m_max = 250
err = (
    []
)  # In the end, it will have the variances of the observables with different block sizes
m_list = []  # The x array for the plot
avg = np.mean(data[N_thermo:, 1:], axis=0)

# Main loop for the Block Average method
while m <= m_max:
    N_blocks = int(N_true / m)
    block_avgs = []  # It will store the averages of each blocks
    ind = 0  # Used for accassing the arrays
    n_bl = 0
    while n_bl < N_blocks:
        aux = np.zeros(4)  # Stores the data before making the average of the block
        for j in range(m):
            aux += data[ind, 1:]
            ind += 1
        if ind > N_thermo:  # The first N_thermo data values are ignored
            block_avgs.append(aux / m)
            n_bl += 1
    block_avgs = np.array(block_avgs)
    total_avg = np.sum(block_avgs, axis=0) / N_blocks
    var2 = np.sum((block_avgs - total_avg) ** 2, axis=0) / ((N_blocks - 1) * N_blocks)
    err.append(np.sqrt(var2))
    m_list.append(m)
    m += 2  # To speed up things a bit. It is unnecessary to go one by one.

# Plot of the results
err = np.array(err)
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(m_list, err[:, 1], lw="1", label="Kinetic", c="#FF474E")
ax.plot(m_list, err[:, 2], lw="1", label="Potential", c="#1982C4")
ax.set_xlim(0)
ax.set_xlabel("Block size")
ax.set_ylabel(r"$\sigma$ (A.U.)")
ax.set_title(r"Variation of $\sigma$ with block size")

plt.legend()
fig.savefig("Var_block.png", format="png", dpi=600)
# plt.show()


# Print the results of the chosen block size after seeing the plot
print("\nThe averages of the results are:")
print("<T'> = {:f}".format(avg[0]))
print("<K'> = {:f}".format(avg[1]))
print("<V'> = {:f}".format(avg[2]))
print("<E'> = {:f}".format(avg[3]))

n_pr = int(
    int(
        input(
            "\nIntroduce the block size to get its results ({} - {}) >> ".format(
                1, m_max
            )
        )
    )
    / 2
)
print("\nVariance calculated with blocks of size {}.".format(m_list[n_pr]))
print("Temperature, Kinetic, Potential, Total")
print(err[n_pr, :])
