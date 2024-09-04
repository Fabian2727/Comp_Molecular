#include <iostream>
#include <vector>
#include <string>
#include <fstream>  // Para leer archivos
#include <sstream>  // Para manejar las líneas del archivo
#include <algorithm>

using namespace std;

int max(int a, int b, int c) {
    return max(max(a, b), c);
}

pair<string, string> needlemanWunsch(const string& seq1, const string& seq2, int match, int mismatch, int gap, int& score) {
    int n = seq1.length();
    int m = seq2.length();

    // Matriz de puntuación
    vector<vector<int>> scoreMatrix(n + 1, vector<int>(m + 1, 0));

    // Inicialización de la matriz
    for (int i = 0; i <= n; i++) {
        scoreMatrix[i][0] = i * gap;
    }
    for (int j = 0; j <= m; j++) {
        scoreMatrix[0][j] = j * gap;
    }

    // Rellenar la matriz
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int scoreDiag = scoreMatrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int scoreLeft = scoreMatrix[i][j - 1] + gap;
            int scoreUp = scoreMatrix[i - 1][j] + gap;
            scoreMatrix[i][j] = max(scoreDiag, scoreLeft, scoreUp);
        }
    }

    // Trazar el alineamiento óptimo
    string alignedSeq1 = "", alignedSeq2 = "";
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch)) {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            i--;
            j--;
        }
        else if (i > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap) {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = "-" + alignedSeq2;
            i--;
        }
        else {
            alignedSeq1 = "-" + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            j--;
        }
    }

    // Asignar el puntaje final
    score = scoreMatrix[n][m];

    // Retornar las secuencias alineadas
    return { alignedSeq1, alignedSeq2 };
}

string cleanSequenceLine(const string& line) {
    string cleaned;
    for (char c : line) {
        if (isalpha(c)) {
            cleaned += c;
        }
    }
    return cleaned;
}

void showSequencesComparisons(const vector<string>& sequences, int match, int mismatch, int gap) {
    for (size_t i = 0; i < sequences.size(); i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            int score;
            auto aligned = needlemanWunsch(sequences[i], sequences[j], match, mismatch, gap, score);
            cout << "[Mejor combinación para la secuencia " << i + 1 << "]: " << aligned.first << endl;
            cout << "[Mejor combinación para la secuencia " << j + 1 << "]: " << aligned.second << endl;
            cout << "Score obtenido: " << score << endl;
            cout << endl;
        }
    }
}

void showScoresOnly(const vector<string>& sequences, int match, int mismatch, int gap) {
    for (size_t i = 0; i < sequences.size(); i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            int score;
            needlemanWunsch(sequences[i], sequences[j], match, mismatch, gap, score);
            cout << "Score para la comparación de la secuencia " << i + 1 << " y la secuencia " << j + 1 << ": " << score << endl;
        }
    }
}

void compareSpecificSequences(int match, int mismatch, int gap) {
    string seq1 = "AAAC";
    string seq2 = "AGC";
    int score;
    auto aligned = needlemanWunsch(seq1, seq2, match, mismatch, gap, score);
    cout << "[Mejor combinación para la secuencia " << seq1 << "]: " << aligned.first << endl;
    cout << "[Mejor combinación para la secuencia " << seq2 << "]: " << aligned.second << endl;
    cout << "Score obtenido: " << score << endl;
}

int main() {
    // Leer secuencias desde el archivo
    ifstream file("C:/Users/Usuario/Desktop/Molecular/Sequencias.txt");
    vector<string> sequences;
    string line;
    string currentSequence = "";

    if (file.is_open()) {
        // Leer el archivo línea por línea
        while (getline(file, line)) {
            // Ignorar líneas vacías y las que comienzan con números (posiciones) o nombres de secuencia
            if (!line.empty() && !isdigit(line[0]) && !isalpha(line[0])) {
                // Limpiar y agregar la secuencia
                currentSequence += cleanSequenceLine(line);
            }
            else if (line.empty()) {
                // Cuando se encuentra una línea vacía, se almacena la secuencia completa
                if (!currentSequence.empty()) {
                    sequences.push_back(currentSequence);
                    currentSequence = "";
                }
            }
        }
        // Agregar la última secuencia al vector si no se agregó
        if (!currentSequence.empty()) {
            sequences.push_back(currentSequence);
        }
        file.close();
    }
    else {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    // Parámetros de puntuación
    int match = 1;
    int mismatch = -1;
    int gap = -2;

    int choice;
    do {
        // Menú
        cout << "Menú de opciones:" << endl;
        cout << "1. Mostrar secuencias comparadas" << endl;
        cout << "2. Mostrar solo puntajes de las comparaciones" << endl;
        cout << "3. Comparar y mostrar el puntaje de AAAC con AGC" << endl;
        cout << "4. Salir" << endl;
        cout << "Elige una opción: ";
        cin >> choice;

        switch (choice) {
        case 1:
            showSequencesComparisons(sequences, match, mismatch, gap);
            break;
        case 2:
            showScoresOnly(sequences, match, mismatch, gap);
            break;
        case 3:
            compareSpecificSequences(match, mismatch, gap);
            break;
        case 4:
            cout << "Saliendo..." << endl;
            break;
        default:
            cout << "Opción no válida. Inténtalo de nuevo." << endl;
            break;
        }

        cout << endl;
    } while (choice != 4);

    return 0;
}
