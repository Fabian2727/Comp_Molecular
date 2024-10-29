import networkx as nx
import matplotlib.pyplot as plt

def alpha(a, b):
    if (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
        return -4
    elif (a == 'C' and b == 'G') or (a == 'G' and b == 'C'):
        return -5
    return 0

def sec_ARN(secuencia):
    L = len(secuencia)
    E = [[0] * L for _ in range(L)]
    P = [[-1] * L for _ in range(L)]

    for d in range(2, L + 1):
        for i in range(L - d + 1):
            j = i + d - 1

            costo = E[i + 1][j - 1] + alpha(secuencia[i], secuencia[j])
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

    for i, j in emparejamientos:
        G.add_edge(i, j)

    cadena_original = [(i, i + 1) for i in range(L - 1)]

    pos = {}
    half = L // 2
    for i in range(half):
        pos[i] = (i, 1)
    for i in range(half, L):
        pos[i] = (L - 1 - i, -1)

    labels = nx.get_node_attributes(G, 'label')

    plt.figure(figsize=(10, 5))
    nx.draw(G, pos, labels=labels, with_labels=True, node_size=700, node_color='lightblue', font_size=10, font_weight='bold')
    nx.draw_networkx_edges(G, pos, edgelist=cadena_original, edge_color='black')
    nx.draw_networkx_edges(G, pos, edgelist=emparejamientos, edge_color='red', style='dashed')

    plt.title('RNA Sequence Graph')
    plt.axis('off')
    plt.show()

def main():
#    secuencia = "GGAAAUCC"
    secuencia = "ACUCGAUUCCGAG"
    E, P = sec_ARN(secuencia)

    emparejamientos = []
    traceback(P, 0, len(secuencia) - 1, secuencia, emparejamientos)

    print("Emparejamientos:", emparejamientos)
    print("Secuencia:", secuencia)
    plot_ARN_grafo(secuencia, emparejamientos)

if __name__ == "__main__":
    main()
