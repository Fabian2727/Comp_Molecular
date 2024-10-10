import matplotlib.pyplot as plt
import networkx as nx

# Función para calcular el solapamiento más largo entre dos secuencias con un linkage t
def overlap(s1, s2, min_length):
    start = 0
    max_overlap = 0
    while True:
        start = s1.find(s2[:min_length], start)
        if start == -1:
            break
        if s1[start:] == s2[:len(s1)-start] and len(s1[start:]) >= min_length:
            max_overlap = len(s1)-start
            break
        start += 1
    return max_overlap

# Función para obtener el complemento reverso de una secuencia de ADN
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_seq = seq[::-1]  # Invertir la secuencia
    reverse_complement_seq = ''.join(complement[base] for base in reversed_seq)
    return reverse_complement_seq

def calculate_consensus_sequence(sequences, target_length=55):
    # Crear una lista para almacenar las secuencias originales y sus complementos reversos
    sequences_with_complement = []
    
    for seq in sequences:
        sequences_with_complement.append(seq)  # Añadir la secuencia original
        reverse_complement_seq = reverse_complement(seq)  # Obtener complemento reverso
        sequences_with_complement.append(reverse_complement_seq)  # Añadir complemento reverso
    
    max_length = max(len(seq) for seq in sequences_with_complement)  # Obtener la longitud máxima de las secuencias
    consensus = ""
    
    for i in range(min(max_length, target_length)):  # Limitar al tamaño objetivo (55 nucleótidos)
        bases = []
        for seq in sequences_with_complement:
            if i < len(seq):  # Solo añadir bases si el índice es válido para la secuencia
                bases.append(seq[i])
        
        if bases:  # Verificar que haya bases disponibles para este índice
            most_common_base = max(set(bases), key=bases.count)
            consensus += most_common_base
        else:
            consensus += "-"  # Si no hay bases, añadir un guion como espacio
            
    return consensus

# Algoritmo greedy para encontrar un camino Hamiltoniano con un valor de linkage t
def hamiltonian_superstring(sequences, t):
    G = nx.DiGraph()
    selected_edges = []
    n = len(sequences)
    in_degree = [0] * n
    out_degree = [0] * n
    disjoint_set = list(range(n))
    
    def find_set(x):
        if disjoint_set[x] != x:
            disjoint_set[x] = find_set(disjoint_set[x])
        return disjoint_set[x]
    
    def union_sets(x, y):
        root_x = find_set(x)
        root_y = find_set(y)
        disjoint_set[root_y] = root_x
    
    # Inicialización
    for i in range(n):
        for j in range(n):
            if i != j:
                overlap_len = overlap(sequences[i], sequences[j], t)
                if overlap_len > 0:
                    G.add_edge(i, j, weight=overlap_len)

    # Ordenar las aristas por peso (mayor solapamiento primero)
    edges = sorted(G.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)

    for edge in edges:
        u, v, weight = edge
        if in_degree[v] == 0 and out_degree[u] == 0 and find_set(u) != find_set(v):
            selected_edges.append((u, v))
            in_degree[v] = 1
            out_degree[u] = 1
            union_sets(u, v)
        if len(selected_edges) == n - 1:  # Si se seleccionan n-1 aristas, termina
            break

    # Visualizar el grafo
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=2000, font_size=10, font_weight='bold', arrows=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title("Grafo de Solapamientos entre Secuencias")
    plt.show()

    # Construir la supercadena final
    final_sequence = sequences[selected_edges[0][0]]
    for u, v in selected_edges:
        overlap_len = overlap(sequences[u], sequences[v], t)
        final_sequence += sequences[v][overlap_len:]
    
    return final_sequence

# Función principal para mostrar el menú
def main_menu():
    
    sequences2 = ["AGTATTGGCAATC", 
                 "CCTTTTGG", 
                 "AATCGATG", 
                 "TTGGCAATCACT", 
                 "ATGCAAACCT"]
    
    sequences1 = [
        "ATCCGTTGAAGCCGCGGGC", 
        "TTAACTCGAGG", 
        "TTAAGTACTGCCCG", 
        "ATCTGTGTCGGG", 
        "CGACTCCCGACACA", 
        "CACAGATCCGTTGAAGCCGCGGG",
        "CTCGAGTTAAGTA", 
        "CGCGGGCAGTACTT"
    ]
    
    while True:
        print("\nMenú:")
        print("1. Secuencia de Consenso")
        print("2. Subgrafos Acíclicos")
        print("3. Salir")
        choice = input("Seleccione una opción (1/2/3): ")
        
        if choice == '1':
                
            print ("USANDO SHORTEST COMMON SUPER STRING: \nInput")

            for seq in sequences1:
                print(seq)

            consensus = calculate_consensus_sequence(sequences1)
            print("\nSecuencia de Consenso:", consensus)
        
        elif choice == '2':
                
            print ("USANDO GREEDY:")
            t = int(input("Ingrese el valor de linkage t: "))
            print("\nSecuencias de Entrada:")
            for seq in sequences2:
                print(seq)
            final_consensus_sequence = hamiltonian_superstring(sequences2, t)  # Cambio aquí
            print("\nSecuencia Consenso Final:", final_consensus_sequence)
        
        elif choice == '3':
            print("Saliendo del programa.")
            break
        
        else:
            print("Opción inválida, por favor intente nuevamente.")

# Ejecutar el menú principal
main_menu()
