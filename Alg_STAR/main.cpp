#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <limits>
#include "needleman.hpp"

using namespace std;

vector<pair<int, pair<string, string>>> StarAlignment(const vector<string>& sequences, vector<vector<int>>& score_matrix) {
    int m = sequences.size();

    if (m == 0) {
        cerr << "Error: No sequences to align." << endl;
        exit(1);
    }

    vector<vector<pair<int, pair<string, string>>>> alignments(m, vector<pair<int, pair<string, string>>>(m, make_pair(0, make_pair(" ", " "))));

    // Comparar todas las secuencias entre sí usando Needleman-Wunsch
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            pair<int, pair<string, string>> alignment = NeedlemanWunsch(sequences[i], sequences[j]);
            alignments[i][j] = alignment;
            alignments[j][i] = alignment;
            score_matrix[i][j] = alignment.first;
            score_matrix[j][i] = alignment.first;
        }
    }

    // Calcular la suma de puntajes para cada secuencia
    vector<int> total_scores(m, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i != j) {
                total_scores[i] += alignments[i][j].first;
            }
        }
    }

    // Encontrar el índice de la secuencia "centro" (la que tiene el mayor puntaje total)
    int max_index = max_element(total_scores.begin(), total_scores.end()) - total_scores.begin();

    // Crear las alineaciones finales con respecto al centro
    vector<pair<int, pair<string, string>>> final_alignments(m, { 0, {" ", " "} });

    for (int i = 0; i < m; i++) {
        if (i != max_index) {
            final_alignments[i] = alignments[i][max_index];
        }
        else {
            final_alignments[i].first = 0;  // El centro no se compara consigo mismo
        }
    }

    return final_alignments;
}

int main() {
    // Leer las secuencias desde el archivo
    vector<pair<string, string>> seqs = read_sequences("BRCA1.txt");

    if (seqs.empty()) {
        cerr << "Error: No sequences were read from the file." << endl;
        return 1;
    }

    // Separar las secuencias forward y reverse
    vector<string> forward_sequences;
    vector<string> reverse_sequences;
    for (const auto& seq : seqs) {
        forward_sequences.push_back(seq.first);   // Forward
        reverse_sequences.push_back(seq.second);  // Reverse
    }

    // Alinear las secuencias Forward
    vector<vector<int>> forward_score_matrix(forward_sequences.size(), vector<int>(forward_sequences.size(), 0));
    vector<pair<int, pair<string, string>>> forward_alignments = StarAlignment(forward_sequences, forward_score_matrix);

    // Alinear las secuencias Reverse
    vector<vector<int>> reverse_score_matrix(reverse_sequences.size(), vector<int>(reverse_sequences.size(), 0));
    vector<pair<int, pair<string, string>>> reverse_alignments = StarAlignment(reverse_sequences, reverse_score_matrix);

    // Imprimir los resultados de las alineaciones Forward
    cout << "Forward Alignments:" << endl;
    for (const auto& alignment : forward_alignments) {
        cout << alignment.second.first << " - " << alignment.second.second << endl;
    }

    // Imprimir la matriz de puntuación Forward
    cout << "Forward Score Matrix:" << endl;
    for (const auto& row : forward_score_matrix) {
        for (int score : row) {
            cout << score << "\t";
        }
        cout << endl;
    }

    // Encontrar la mejor secuencia forward
    int best_forward_index = max_element(forward_alignments.begin(), forward_alignments.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; }) - forward_alignments.begin();

    cout << "Best Forward Sequence:" << endl;
    cout << forward_sequences[best_forward_index] << endl;

    // Comparar la mejor secuencia forward con todas las demás
    cout << "Best Forward Sequence vs All Others:" << endl;
    for (int i = 0; i < forward_sequences.size(); i++) {
        if (i != best_forward_index) {
            auto alignment = NeedlemanWunsch(forward_sequences[best_forward_index], forward_sequences[i]);
            cout << forward_sequences[best_forward_index] << " vs " << forward_sequences[i] << ":" << endl;
            cout << alignment.second.first << " - " << alignment.second.second << endl;
        }
    }

    // Imprimir los resultados de las alineaciones Reverse
    cout << "Reverse Alignments:" << endl;
    for (const auto& alignment : reverse_alignments) {
        cout << alignment.second.first << " - " << alignment.second.second << endl;
    }

    // Imprimir la matriz de puntuación Reverse
    cout << "Reverse Score Matrix:" << endl;
    for (const auto& row : reverse_score_matrix) {
        for (int score : row) {
            cout << score << "\t";
        }
        cout << endl;
    }

    // Encontrar la mejor secuencia reverse
    int best_reverse_index = max_element(reverse_alignments.begin(), reverse_alignments.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; }) - reverse_alignments.begin();

    cout << "Best Reverse Sequence:" << endl;
    cout << reverse_sequences[best_reverse_index] << endl;

    // Comparar la mejor secuencia reverse con todas las demás
    cout << "Best Reverse Sequence vs All Others:" << endl;
    for (int i = 0; i < reverse_sequences.size(); i++) {
        if (i != best_reverse_index) {
            auto alignment = NeedlemanWunsch(reverse_sequences[best_reverse_index], reverse_sequences[i]);
            cout << reverse_sequences[best_reverse_index] << " vs " << reverse_sequences[i] << ":" << endl;
            cout << alignment.second.first << " - " << alignment.second.second << endl;
        }
    }

    return 0;
}
