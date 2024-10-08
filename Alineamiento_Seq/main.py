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
    reverse_seq = seq[::-1]
    return ''.join([complement[base] for base in reverse_seq])

# Algoritmo para encontrar la secuencia de consenso
def calculate_consensus_sequence(sequences):
    max_length = max(len(seq) for seq in sequences)  # Obtener la longitud máxima de las secuencias
    consensus = ""
    
    for i in range(max_length):
        bases = []
        for seq in sequences:
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
    while len(sequences) > 1:
        max_overlap_len = -1
        best_pair = (0, 0)
        best_merged_seq = ""
        
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i != j:
                    overlap_len = overlap(sequences[i], sequences[j], t)
                    if overlap_len > max_overlap_len:
                        max_overlap_len = overlap_len
                        best_pair = (i, j)
                        best_merged_seq = sequences[i] + sequences[j][overlap_len:]
                    
                    rev_seq = reverse_complement(sequences[j])
                    overlap_len_rev = overlap(sequences[i], rev_seq, t)
                    if overlap_len_rev > max_overlap_len:
                        max_overlap_len = overlap_len_rev
                        best_pair = (i, j)
                        best_merged_seq = sequences[i] + rev_seq[overlap_len_rev:]
        
        if max_overlap_len < t:
            print("No se puede continuar con el ensamblaje. No hay más solapamientos con al menos", t, "bases.")
            break
        
        i, j = best_pair
        sequences[i] = best_merged_seq
        sequences.pop(j)
        
        G.add_edge(f'seq{i}', f'seq{j}', weight=max_overlap_len)
    
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=2000, font_size=10, font_weight='bold', arrows=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title("Grafo de Solapamientos entre Secuencias")
    plt.show()
    
    return sequences[0]

# Función principal para mostrar el menú
def main_menu():
    sequences = [
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
            consensus = calculate_consensus_sequence(sequences)  # Cambio aquí
            print("\nSecuencia de Consenso:", consensus)
        
        elif choice == '2':
            t = int(input("Ingrese el valor de linkage t: "))
            print("\nSecuencias de Entrada:")
            for seq in sequences:
                print(seq)
            final_consensus_sequence = hamiltonian_superstring(sequences, t)  # Cambio aquí
            print("\nSecuencia Consenso Final:", final_consensus_sequence)
        
        elif choice == '3':
            print("Saliendo del programa.")
            break
        
        else:
            print("Opción inválida, por favor intente nuevamente.")

# Ejecutar el menú principal
main_menu()
