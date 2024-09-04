import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Función para crear la matriz de coincidencias
def dot_matrix(string1, string2):
    matrix = np.zeros((len(string1), len(string2)))

    for i in range(len(string1)):
        for j in range(len(string2)):
            if string1[i] == string2[j]:
                matrix[i][j] = 1

    return matrix

# Función para mostrar la matriz con un mapa de colores en blanco y negro
def plot_dot_matrix(matrix, label_x, label_y):
    # Definir un mapa de colores personalizado: blanco para no coincidencias (0), negro para coincidencias (1)
    cmap = ListedColormap(['white', 'black'])

    plt.imshow(matrix, cmap=cmap, interpolation='nearest')

    # Configurar etiquetas en los ejes X e Y, sin las secuencias, solo los nombres
    plt.xticks([])  # Eliminar las etiquetas de las secuencias en el eje X
    plt.yticks([])  # Eliminar las etiquetas de las secuencias en el eje Y
    plt.xlabel(label_x)  # Etiqueta del eje X
    plt.ylabel(label_y)  # Etiqueta del eje Y
    plt.title('Dot Matrix of Strings')
    plt.show()

# Leer las secuencias desde el archivo
sequences = []
with open('C:/Users/Usuario/Desktop/Molecular/newSeq.txt', 'r') as file:
    for line in file:
        sequence = line.strip()
        sequences.append(sequence)

bacteria = sequences[0]
sarscov = sequences[1]
influenza = sequences[2]

# Crear y mostrar la matriz de coincidencias
matrix = dot_matrix(sarscov, bacteria)
plot_dot_matrix(matrix, 'Influenza', 'SarsCov')
