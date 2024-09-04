#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

// Función para calcular el máximo de tres números
int maximo(int a, int b, int c) {
    return max(max(a, b), c);
}

// Función para contar el número de rupturas en una secuencia alineada
int contarRupturas(const string& secuenciaAlineada) {
    int rupturas = 0;
    bool enBloqueGuiones = false;

    for (char c : secuenciaAlineada) {
        if (c == '-') {
            if (!enBloqueGuiones) {
                rupturas++;
                enBloqueGuiones = true;
            }
        }
        else {
            enBloqueGuiones = false;
        }
    }

    return rupturas;
}

// Función para alinear dos secuencias usando Needleman-Wunsch
pair<string, string> alineacionNeedlemanWunsch(const string& secuencia1, const string& secuencia2, int coincidencia, int desajuste, int gap, int& puntaje) {
    int n = secuencia1.length();
    int m = secuencia2.length();

    // Matriz de puntuación
    vector<vector<int>> matrizPuntaje(n + 1, vector<int>(m + 1, 0));

    // Inicialización de la matriz
    for (int i = 0; i <= n; i++) {
        matrizPuntaje[i][0] = i * gap;
    }
    for (int j = 0; j <= m; j++) {
        matrizPuntaje[0][j] = j * gap;
    }

    // Rellenar la matriz
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int puntajeDiagonal = matrizPuntaje[i - 1][j - 1] + (secuencia1[i - 1] == secuencia2[j - 1] ? coincidencia : desajuste);
            int puntajeIzquierda = matrizPuntaje[i][j - 1] + gap;
            int puntajeArriba = matrizPuntaje[i - 1][j] + gap;
            matrizPuntaje[i][j] = maximo(puntajeDiagonal, puntajeIzquierda, puntajeArriba);
        }
    }

    // Trazar el alineamiento óptimo
    string secuenciaAlineada1 = "", secuenciaAlineada2 = "";
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && matrizPuntaje[i][j] == matrizPuntaje[i - 1][j - 1] + (secuencia1[i - 1] == secuencia2[j - 1] ? coincidencia : desajuste)) {
            secuenciaAlineada1 = secuencia1[i - 1] + secuenciaAlineada1;
            secuenciaAlineada2 = secuencia2[j - 1] + secuenciaAlineada2;
            i--;
            j--;
        }
        else if (i > 0 && matrizPuntaje[i][j] == matrizPuntaje[i - 1][j] + gap) {
            secuenciaAlineada1 = secuencia1[i - 1] + secuenciaAlineada1;
            secuenciaAlineada2 = "-" + secuenciaAlineada2;
            i--;
        }
        else {
            secuenciaAlineada1 = "-" + secuenciaAlineada1;
            secuenciaAlineada2 = secuencia2[j - 1] + secuenciaAlineada2;
            j--;
        }
    }

    // Asignar el puntaje final
    puntaje = matrizPuntaje[n][m];

    // Retornar las secuencias alineadas
    return { secuenciaAlineada1, secuenciaAlineada2 };
}

// Función para limpiar una línea de secuencia
string limpiarSecuencia(const string& linea) {
    string secuenciaLimpia;
    for (char c : linea) {
        if (isalpha(c)) {
            secuenciaLimpia += c;
        }
    }
    return secuenciaLimpia;
}

// Función para mostrar el alineamiento con menos rupturas
void mostrarAlineacionConMenosRupturas(const vector<string>& secuencias, int coincidencia, int desajuste, int gap) {
    if (secuencias.size() < 3) {
        cout << "Se necesitan al menos 3 secuencias para realizar esta operación." << endl;
        return;
    }

    // Comparar secuencia 1 con secuencia 2
    int puntaje1_2;
    auto alineacion1_2 = alineacionNeedlemanWunsch(secuencias[0], secuencias[1], coincidencia, desajuste, gap, puntaje1_2);
    int rupturas1_2 = contarRupturas(alineacion1_2.first) + contarRupturas(alineacion1_2.second);

    // Comparar secuencia 1 con secuencia 3
    int puntaje1_3;
    auto alineacion1_3 = alineacionNeedlemanWunsch(secuencias[0], secuencias[2], coincidencia, desajuste, gap, puntaje1_3);
    int rupturas1_3 = contarRupturas(alineacion1_3.first) + contarRupturas(alineacion1_3.second);

    // Comparar secuencia 2 con secuencia 3
    int puntaje2_3;
    auto alineacion2_3 = alineacionNeedlemanWunsch(secuencias[1], secuencias[2], coincidencia, desajuste, gap, puntaje2_3);
    int rupturas2_3 = contarRupturas(alineacion2_3.first) + contarRupturas(alineacion2_3.second);

    // Mostrar los resultados
    cout << "Comparación de secuencia 1 con secuencia 2:" << endl;
    cout << alineacion1_2.first << endl;
    cout << alineacion1_2.second << endl;
    cout << "Número de rupturas: " << rupturas1_2 << endl << endl;

    cout << "Comparación de secuencia 1 con secuencia 3:" << endl;
    cout << alineacion1_3.first << endl;
    cout << alineacion1_3.second << endl;
    cout << "Número de rupturas: " << rupturas1_3 << endl << endl;

    cout << "Comparación de secuencia 2 con secuencia 3:" << endl;
    cout << alineacion2_3.first << endl;
    cout << alineacion2_3.second << endl;
    cout << "Número de rupturas: " << rupturas2_3 << endl;
}

// Función para mostrar las comparaciones de secuencias
void mostrarComparacionesSecuencias(const vector<string>& secuencias, int coincidencia, int desajuste, int gap) {
    for (size_t i = 0; i < secuencias.size(); i++) {
        for (size_t j = i + 1; j < secuencias.size(); j++) {
            int puntaje;
            auto alineacion = alineacionNeedlemanWunsch(secuencias[i], secuencias[j], coincidencia, desajuste, gap, puntaje);
            cout << "[Mejor combinación para la secuencia " << i + 1 << "]: " << alineacion.first << endl;
            cout << "[Mejor combinación para la secuencia " << j + 1 << "]: " << alineacion.second << endl;
            cout << "Puntaje obtenido: " << puntaje << endl;
            cout << endl;
        }
    }
}

// Función para mostrar solo los puntajes de las comparaciones
void mostrarSoloPuntajes(const vector<string>& secuencias, int coincidencia, int desajuste, int gap) {
    for (size_t i = 0; i < secuencias.size(); i++) {
        for (size_t j = i + 1; j < secuencias.size(); j++) {
            int puntaje;
            alineacionNeedlemanWunsch(secuencias[i], secuencias[j], coincidencia, desajuste, gap, puntaje);
            cout << "Puntaje para la comparación de la secuencia " << i + 1 << " y la secuencia " << j + 1 << ": " << puntaje << endl;
        }
    }
}

// Función para comparar secuencias específicas (AAAC y AGC)
void compararSecuenciasEspecificas(int coincidencia, int desajuste, int gap) {
    string secuencia1 = "AAAC";
    string secuencia2 = "AGC";
    int puntaje;
    auto alineacion = alineacionNeedlemanWunsch(secuencia1, secuencia2, coincidencia, desajuste, gap, puntaje);
    cout << "[Mejor combinación para la secuencia " << secuencia1 << "]: " << alineacion.first << endl;
    cout << "[Mejor combinación para la secuencia " << secuencia2 << "]: " << alineacion.second << endl;
    cout << "Puntaje obtenido: " << puntaje << endl;
}

int main() {
    // Leer secuencias desde el archivo
    ifstream archivo("C:/Users/Usuario/Desktop/Molecular/Sequencias.txt");
    vector<string> secuencias;
    string linea;
    string secuenciaActual = "";

    if (archivo.is_open()) {
        // Leer el archivo línea por línea
        while (getline(archivo, linea)) {
            // Ignorar líneas vacías y las que comienzan con números o nombres de secuencia
            if (!linea.empty() && !isdigit(linea[0]) && !isalpha(linea[0])) {
                secuenciaActual += limpiarSecuencia(linea);
            }
            else if (linea.empty()) {
                if (!secuenciaActual.empty()) {
                    secuencias.push_back(secuenciaActual);
                    secuenciaActual = "";
                }
            }
        }
        if (!secuenciaActual.empty()) {
            secuencias.push_back(secuenciaActual);
        }
        archivo.close();
    }
    else {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    // Parámetros de puntuación
    int coincidencia = 1;
    int desajuste = -1;
    int gap = -2;

    int opcion;
    do {
        // Menú
        cout << "Menú de opciones:" << endl;
        cout << "1. Mostrar secuencias comparadas" << endl;
        cout << "2. Mostrar solo puntajes de las comparaciones" << endl;
        cout << "3. Comparar y mostrar el puntaje de AAAC con AGC" << endl;
        cout << "4. Mostrar la alineación con menos rupturas" << endl;
        cout << "5. Salir" << endl;
        cout << "Elige una opción: ";
        cin >> opcion;

        switch (opcion) {
        case 1:
            mostrarComparacionesSecuencias(secuencias, coincidencia, desajuste, gap);
            break;
        case 2:
            mostrarSoloPuntajes(secuencias, coincidencia, desajuste, gap);
            break;
        case 3:
            compararSecuenciasEspecificas(coincidencia, desajuste, gap);
            break;
        case 4:
            mostrarAlineacionConMenosRupturas(secuencias, coincidencia, desajuste, gap);
            break;
        case 5:
            cout << "Saliendo..." << endl;
            break;
        default:
            cout << "Opción no válida. Inténtalo de nuevo." << endl;
            break;
        }

        cout << endl;
    } while (opcion != 5);

    return 0;
}
