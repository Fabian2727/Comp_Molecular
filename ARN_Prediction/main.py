import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

def alpha(a, b, energy_function=1):
    if energy_function == 1:
        if (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
            return -4
        elif (a == 'C' and b == 'G') or (a == 'G' and b == 'C'):
            return -5
        elif (a == 'G' and b == 'U') or (a == 'U' and b == 'G'):
            return -1
    elif energy_function == 2:
        if (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
            return -1
        elif (a == 'C' and b == 'G') or (a == 'G' and b == 'C'):
            return -1
    return 0

def is_complementary(a, b):
    return (a == 'A' and b == 'U') or (a == 'U' and b == 'A') or \
           (a == 'C' and b == 'G') or (a == 'G' and b == 'C') or \
           (a == 'G' and b == 'U') or (a == 'U' and b == 'G')

def sec_ARN(secuencia, energy_function=1):
    L = len(secuencia)
    E = [[0] * L for _ in range(L)]
    P = [[-1] * L for _ in range(L)]

    for d in range(2, L + 1):
        for i in range(L - d + 1):
            j = i + d - 1

            costo = E[i + 1][j - 1] + alpha(secuencia[i], secuencia[j], energy_function)
            if costo < E[i][j - 1]:
                E[i][j] = costo
                P[i][j] = -2
            else:
                E[i][j] = E[i][j - 1]
                P[i][j] = -1

            for k in range(i + 1, j):
                if E[i][j] > E[i][k] + E[k + 1][j]:
                    E[i][j] = E[i][k] + E[k + 1][j]
                    P[i][j] = k

    return E, P

def traceback(P, i, j, secuencia, emparejamientos):
    if i >= j:
        return
    if P[i][j] == -1:
        traceback(P, i, j - 1, secuencia, emparejamientos)
    elif P[i][j] == -2:
        if j != i + 1:
            emparejamientos.append((i, j))
        traceback(P, i + 1, j - 1, secuencia, emparejamientos)
    else:
        k = P[i][j]
        traceback(P, i, k, secuencia, emparejamientos)
        traceback(P, k + 1, j, secuencia, emparejamientos)

def plot_ARN_grafo(secuencia, emparejamientos):
    G = nx.Graph()
    L = len(secuencia)

    for i in range(L):
        G.add_node(i, label=secuencia[i])

    # Separar los emparejamientos complementarios para colorearlos en verde
    complementarios = [(i, j) for i, j in emparejamientos if is_complementary(secuencia[i], secuencia[j])]
    no_complementarios = [pair for pair in emparejamientos if pair not in complementarios]

    # Enlaces de la cadena original
    cadena_original = [(i, i + 1) for i in range(L - 1)]

    # Distribuir nodos en un círculo
    pos = {i: (np.cos(2 * np.pi * i / L), np.sin(2 * np.pi * i / L)) for i in range(L)}

    labels = nx.get_node_attributes(G, 'label')

    plt.figure(figsize=(10, 10))
    nx.draw(G, pos, labels=labels, with_labels=True, node_size=700, node_color='lightblue', font_size=10, font_weight='bold')
    
    # Dibujar los enlaces de la cadena en negro
    nx.draw_networkx_edges(G, pos, edgelist=cadena_original, edge_color='black')
    
    # Dibujar los emparejamientos no complementarios en rojo
    nx.draw_networkx_edges(G, pos, edgelist=no_complementarios, edge_color='red', style='dashed')
    
    # Dibujar los emparejamientos complementarios en verde
    nx.draw_networkx_edges(G, pos, edgelist=complementarios, edge_color='green', style='dashed')

    plt.title('RNA Sequence Graph with Complementary Pairings in Green')
    plt.axis('off')
    plt.show()

def main():
    print("Select energy function:")
    print("1. Function 1")
    print("2. Function 2")
    energy_choice = int(input("Enter your choice (1 or 2): "))
    
    print("Select RNA sequence:")
    print("1. GGAAAUCC")
    print("2. ACUCGAUUCCGAG")
    sequence_choice = int(input("Enter your choice (1 or 2): "))
    
    secuencia = "GGAAAUCC" if sequence_choice == 1 else "ACUCGAUUCCGAG"
    
    E, P = sec_ARN(secuencia, energy_function=energy_choice)

    emparejamientos = []
    traceback(P, 0, len(secuencia) - 1, secuencia, emparejamientos)

    print("Emparejamientos:", emparejamientos)
    print("Secuencia:", secuencia)
    plot_ARN_grafo(secuencia, emparejamientos)

if __name__ == "__main__":
    main()
