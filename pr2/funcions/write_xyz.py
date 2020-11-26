def write_file(particles, filetype):
    f = open(filetype, "a")
    f.write(str(len(particles)) + "\n\n")
    for particle in particles:
        string = ""
        for j in particle:
            string += str(j) + " "
        f.write("C " + string + "\n")
    f.close
