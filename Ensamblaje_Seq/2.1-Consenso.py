def reverse_complement(s):
    complements = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return ''.join(complements[nucleotide] for nucleotide in reversed(s))

def find_overlap(seq1, seq2):
    overlap = 0
    length = min(len(seq1), len(seq2))
    for i in range(length):
        if seq1[-length + i:] == seq2[:length - i]:
            overlap = length - i
            break
    return overlap

def greedy_hamiltonian_path(sequences):
    visited = [False] * len(sequences)
    hamiltonian_path = []
    
    # Create overlap matrix
    overlap_matrix = [[0 for _ in range(len(sequences))] for _ in range(len(sequences))]
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i != j:
                overlap_matrix[i][j] = find_overlap(sequences[i], sequences[j])

    # Find the sequence with the maximum overlap to start
    max_overlap = 0
    start_vertex = -1
    for i in range(len(sequences)):
        current_overlap = 0
        for j in range(len(sequences)):
            if i != j and overlap_matrix[i][j] > current_overlap:
                current_overlap = overlap_matrix[i][j]
        if current_overlap > max_overlap:
            max_overlap = current_overlap
            start_vertex = i

    if start_vertex == -1:
        start_vertex = 0

    hamiltonian_path.append((0, sequences[start_vertex]))
    visited[start_vertex] = True

    # Build the Hamiltonian path
    for _ in range(len(sequences) - 1):
        max_overlap = 0
        next_vertex = -1
        for j in range(len(sequences)):
            if not visited[j]:
                current_sequence = sequences[j]
                rev_complement = reverse_complement(current_sequence)
                overlap = max(find_overlap(hamiltonian_path[-1][1], current_sequence), 
                              find_overlap(hamiltonian_path[-1][1], rev_complement))
                if overlap > max_overlap:
                    max_overlap = overlap
                    next_vertex = j
        
        if next_vertex != -1:
            next_sequence = sequences[next_vertex]
            rev_complement = reverse_complement(next_sequence)
            if find_overlap(hamiltonian_path[-1][1], next_sequence) >= find_overlap(hamiltonian_path[-1][1], rev_complement):
                hamiltonian_path.append((max_overlap, next_sequence))
            else:
                hamiltonian_path.append((max_overlap, rev_complement))
            visited[next_vertex] = True

    return hamiltonian_path

def main():
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

    hamiltonian_path = greedy_hamiltonian_path(sequences1)

    alignment = ""
    consensus = hamiltonian_path[0][1]  # Iniciar la secuencia consenso con la primera secuencia completa
    s_ant = 0
    pad = 0

    for i, (overlap, sequence) in enumerate(hamiltonian_path):
        if i > 0:
            # Solo añadir la parte que no se solapa a la secuencia consenso
            consensus += sequence[overlap:]
        
        # Alineamiento para impresión
        pad += s_ant - overlap
        alignment += " " * pad + sequence + '\n'
        s_ant = len(sequence)

    print("Camino Hamiltoniano:")
    for overlap, sequence in hamiltonian_path:
        print(f"Overlap: {overlap}, Secuencia visitada: {sequence}")

    print("\nAlineamiento:")
    print(alignment)

    # Imprimir la secuencia consenso final
    print("-----------------------------------------")
    print(consensus)

if __name__ == "__main__":
    main()
