import numpy as np
import math

"""sequences = {
    'S1': "ATTGCCATT",
    'S2': "ATGGCCATT",
    'S3': "ATCCAATTTT",
    'S4': "ATCTTCTT",
    'S5': "ACTGACC"
}

"""
sequences = {
    'S1': "CONCHA",
    'S2': "FLORES",
    'S3': "PINO",
    'S4': "QUISPE",
    'S5': "SALAS",
    'S6': "TUPAC",
    'S7': "VILLANUEVA"
}

# Función para calcular la distancia de Hamming entre dos secuencias
def hamming_distance(seq1, seq2):
    length = min(len(seq1), len(seq2))
    differences = sum(1 for i in range(length) if seq1[i] != seq2[i])
    return differences, length

# Función para calcular la matriz de distancias corregida usando Ecuación 1
def calculate_distance_matrix(sequences):
    keys = list(sequences.keys())
    n = len(keys)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            diffs, length = hamming_distance(sequences[keys[i]], sequences[keys[j]])
            p = diffs / length
            # Aplicamos la ecuación de tasa de sustitución
            if p < 0.75:
                dist = -0.75 * math.log(1 - (4/3) * p)
            else:
                dist = 1  # En lugar de infinito, asignamos una distancia máxima de 1
            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist

    return dist_matrix, keys

# Calculamos la matriz de distancias
dist_matrix, seq_keys = calculate_distance_matrix(sequences)

# Imprime la matriz de distancias
print("Matriz de distancias corregida:")
print(dist_matrix)

# Algoritmo UPGMA para construcción de árbol filogenético
def upgma(dist_matrix, keys):
    clusters = {i: [key] for i, key in enumerate(keys)}
    heights = {i: 0 for i in range(len(keys))}  # Alturas para la unión de nodos
    
    while len(clusters) > 1:
        # Encontrar los clusters más cercanos
        min_dist = float('inf')
        to_merge = None
        for i in clusters:
            for j in clusters:
                if i < j and dist_matrix[i][j] < min_dist:
                    min_dist = dist_matrix[i][j]
                    to_merge = (i, j)
        
        # Verificar si encontramos un par de clusters a unir
        if to_merge is None:
            break
        
        # Unir dos clusters
        i, j = to_merge
        new_cluster = clusters[i] + clusters[j]
        print(f"Uniendo {clusters[i]} y {clusters[j]} con distancia {min_dist}")
        del clusters[j]
        clusters[i] = new_cluster

        # Actualizar las alturas para las uniones
        heights[i] = min_dist / 2

        # Actualizar la matriz de distancias
        for k in clusters:
            if k != i:
                dist_matrix[i][k] = dist_matrix[k][i] = (dist_matrix[i][k] + dist_matrix[j][k]) / 2
        dist_matrix[i][i] = 0
    
    return clusters
    
upgma_tree = upgma(dist_matrix, seq_keys)
print("Árbol filogenético UPGMA:", upgma_tree)

# Algoritmo Neighbor Joining para construcción de árbol filogenético
def neighbor_joining(dist_matrix, keys):
    clusters = {i: [key] for i, key in enumerate(keys)}
    
    while len(clusters) > 2:
        n = len(clusters)
        total_distances = np.sum(dist_matrix, axis=1)

        # Calcula la corrección Q
        q_matrix = np.zeros_like(dist_matrix)
        for i in clusters:
            for j in clusters:
                if i != j:
                    q_matrix[i][j] = (n - 2) * dist_matrix[i][j] - total_distances[i] - total_distances[j]

        # Encontrar los pares más cercanos
        i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)

        # Unir los dos clusters más cercanos
        new_cluster = clusters[i] + clusters[j]
        print(f"Uniendo {clusters[i]} y {clusters[j]} con distancia {dist_matrix[i][j]}")
        del clusters[j]
        clusters[i] = new_cluster

        # Actualizar la matriz de distancias
        for k in clusters:
            if k != i:
                dist_matrix[i][k] = dist_matrix[k][i] = (dist_matrix[i][k] + dist_matrix[j][k] - dist_matrix[i][j]) / 2
        dist_matrix[i][i] = 0
    
    return clusters

nj_tree = neighbor_joining(dist_matrix, seq_keys)
print("Árbol filogenético Neighbor Joining:", nj_tree)
