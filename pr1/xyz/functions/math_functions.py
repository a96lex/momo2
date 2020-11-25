def vector_module(array):
    suma = 0
    for i in array:
        suma += i
    return suma


def suma_vectors(array1, array2):
    suma = []
    for i in range(len(array1)):
        suma.append(array1[i] + array2[i])
    return suma


def suma_vectors_prod(array1, s1, array2, s2):
    suma = []
    for i in range(len(array1)):
        suma.append(array1[i] * s1 + array2[i] * s2)
    return suma


def scalar_prod_vectors(array1, scalar):
    suma = []
    for i in range(len(array1)):
        suma.append(array1[i] * scalar)
    return suma
