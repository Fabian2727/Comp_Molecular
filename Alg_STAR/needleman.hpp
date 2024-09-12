#include <iostream>
#include <string>
#include <stack>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cctype>
#include <algorithm>
using namespace std;

/* Needleman-Wunsch */
pair<int, vector<vector<vector<short>>>> NeedlemanWunsch_(string s1, string s2) {
    if (s1.empty() || s2.empty()) {
        cerr << "Error: One of the input sequences is empty." << endl;
        exit(1);
    }

    vector<vector<vector<short>>> paths;
    vector<vector<int>> score;
    int m = s1.length();
    int n = s2.length();

    score.assign(m + 1, vector<int>(n + 1, 0));
    paths.assign(m + 1, vector<vector<short>>(n + 1, vector<short>()));

    score[0][0] = 0;

    for (int i = 1; i <= m; i++) {
        score[i][0] = score[i - 1][0] - 2;
        paths[i][0].push_back(1);
    }

    for (int j = 1; j <= n; j++) {
        score[0][j] = score[0][j - 1] - 2;
        paths[0][j].push_back(2);
    }

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int score_a = score[i - 1][j - 1] + (s1[i - 1] == s2[j - 1]) * 2 - 1;
            int score_b = score[i - 1][j] - 2;
            int score_c = score[i][j - 1] - 2;
            int max_score = max({ score_a, score_b, score_c });

            if (score_a == max_score)
                paths[i][j].push_back(0);
            if (score_b == max_score)
                paths[i][j].push_back(1);
            if (score_c == max_score)
                paths[i][j].push_back(2);

            score[i][j] = max_score;
        }
    }

    return { score[m][n], paths };
}

pair<string, string> getBestAlignment(int i, int j, const string& s1, const string& s2, const vector<vector<vector<short>>>& paths) {
    if (i == 0 && j == 0)
        return { "", "" };

    pair<string, string> alignment;

    if (paths[i][j].empty()) {
        cerr << "Error: No valid path at position (" << i << ", " << j << ")" << endl;
        exit(1);
    }

    short path = paths[i][j][0];

    if (path == 0) {
        alignment = getBestAlignment(i - 1, j - 1, s1, s2, paths);
        alignment.first += s1[i - 1];
        alignment.second += s2[j - 1];
    }
    else if (path == 1) {
        alignment = getBestAlignment(i - 1, j, s1, s2, paths);
        alignment.first += s1[i - 1];
        alignment.second += '-';
    }
    else if (path == 2) {
        alignment = getBestAlignment(i, j - 1, s1, s2, paths);
        alignment.first += '-';
        alignment.second += s2[j - 1];
    }
    return alignment;
}

pair<int, pair<string, string>> NeedlemanWunsch(string s1, string s2) {
    if (s1.empty() || s2.empty()) {
        cerr << "Error: Empty sequence provided to Needleman-Wunsch." << endl;
        exit(1);
    }

    pair<int, vector<vector<vector<short>>>> s_a = NeedlemanWunsch_(s1, s2);
    int m = s1.length(), n = s2.length();
    pair<string, string> alignment = getBestAlignment(m, n, s1, s2, s_a.second);
    return { s_a.first, alignment };
}

vector<pair<string, string>> read_sequences(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }

    vector<pair<string, string>> sequences;
    string line;
    string forward, reverse;

    while (getline(file, line)) {
        if (line.empty()) continue;

        // Identificar líneas de secuencias forward y reverse
        if (line.find("F:") != string::npos) {
            forward = line.substr(line.find("5-") + 2, line.find("-3") - line.find("5-") - 2);
        }
        else if (line.find("R:") != string::npos) {
            reverse = line.substr(line.find("5-") + 2, line.find("-3") - line.find("5-") - 2);
            // Asumir que después de leer un reverse, se ha leído un forward correspondiente
            if (!forward.empty() && !reverse.empty()) {
                sequences.push_back({ forward, reverse });
                forward.clear(); // Limpiar para la próxima secuencia
                reverse.clear(); // Limpiar para la próxima secuencia
            }
        }
    }

    file.close();

    if (sequences.size() != 6) {
        cerr << "Error: Expected 6 sequences, but found " << sequences.size() << endl;
        exit(1);
    }

    return sequences;
}
