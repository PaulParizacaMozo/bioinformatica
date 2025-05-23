#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Parámetros para el alineamiento
const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

// Estructura para almacenar la información detallada de un alineamiento local
struct AlineamientoInfo {
  string s1_alineada;
  string s2_alineada;
  int start_s1;
  int end_s1;
  int start_s2;
  int end_s2;
};

// Estructura para los resultados del alineamiento local
struct ResultadoAlineamientoLocal {
  int scoreMayor;
  vector<vector<int>> matrizScores;
  vector<AlineamientoInfo> alineamientos;
};

// Función para reconstruir un alineamiento local
void reconstruir(const string &s1, const string &s2, const vector<vector<int>> &matriz, int end_row, int end_col,
                 vector<AlineamientoInfo> &todosLosAlineamientos) {

  if (end_row == 0 || end_col == 0 || matriz[end_row][end_col] == 0) {
    return;
  }

  string alin1_rev, alin2_rev;
  int i = end_row;
  int j = end_col;

  // Reconstruir el alineamiento hacia atrás
  while (i > 0 && j > 0 && matriz[i][j] != 0) {
    int scoreActual = matriz[i][j];
    // El carácter actual de s1 es s1[i-1], de s2 es s2[j-1]
    int scoreDiagonal = (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH);

    if (scoreActual == matriz[i - 1][j - 1] + scoreDiagonal) {
      alin1_rev += s1[i - 1];
      alin2_rev += s2[j - 1];
      i--;
      j--;
    } else if (scoreActual == matriz[i - 1][j] + GAP) {
      alin1_rev += s1[i - 1];
      alin2_rev += '-';
      i--;
    } else if (scoreActual == matriz[i][j - 1] + GAP) {
      alin1_rev += '-';
      alin2_rev += s2[j - 1];
      j--;
    } else {
      break;
    }
  }

  AlineamientoInfo info;
  reverse(alin1_rev.begin(), alin1_rev.end());
  reverse(alin2_rev.begin(), alin2_rev.end());
  info.s1_alineada = alin1_rev;
  info.s2_alineada = alin2_rev;

  info.end_s1 = end_row - 1;
  info.end_s2 = end_col - 1;
  info.start_s1 = i;
  info.start_s2 = j;

  if (!info.s1_alineada.empty() || !info.s2_alineada.empty()) {
    bool duplicado = false;
    for (const auto &existente : todosLosAlineamientos) {
      if (existente.s1_alineada == info.s1_alineada && existente.s2_alineada == info.s2_alineada &&
          existente.start_s1 == info.start_s1 && existente.end_s1 == info.end_s1 && existente.start_s2 == info.start_s2 &&
          existente.end_s2 == info.end_s2) {
        duplicado = true;
        break;
      }
    }
    if (!duplicado) {
      todosLosAlineamientos.push_back(info);
    }
  }
}

// Implementación del alineamiento local
ResultadoAlineamientoLocal alineamientoLocal(const string &s1, const string &s2) {
  int n = s1.length();
  int m = s2.length();

  vector<vector<int>> matriz(n + 1, vector<int>(m + 1));
  int scoreMayor = 0;
  vector<pair<int, int>> celdasMaxScore; // Almacena coordenadas de celdas con scoreMayor

  for (int i = 0; i <= n; ++i)
    matriz[i][0] = 0;
  for (int j = 0; j <= m; ++j)
    matriz[0][j] = 0;

  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int scoreDiagonal = matriz[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH);
      int scoreArriba = matriz[i - 1][j] + GAP;
      int scoreIzquierda = matriz[i][j - 1] + GAP;
      matriz[i][j] = max({0, scoreDiagonal, scoreArriba, scoreIzquierda});

      if (matriz[i][j] > scoreMayor) {
        scoreMayor = matriz[i][j];
        celdasMaxScore.clear();
        celdasMaxScore.push_back({i, j});
      } else if (matriz[i][j] == scoreMayor && scoreMayor > 0) {
        celdasMaxScore.push_back({i, j});
      }
    }
  }

  ResultadoAlineamientoLocal resultado;
  resultado.scoreMayor = scoreMayor;
  resultado.matrizScores = matriz;

  // reconstruccion con el score mayor
  for (const auto &celda : celdasMaxScore) {
    reconstruir(s1, s2, matriz, celda.first, celda.second, resultado.alineamientos);
  }
  return resultado;
}

// Función guardar resultados
void guardarResultados(const string &nombreArchivo, const ResultadoAlineamientoLocal &resultado) {
  ofstream archivoSalida(nombreArchivo);
  if (archivoSalida.is_open()) {
    archivoSalida << "* Score final(Optimo): " << resultado.scoreMayor << endl;

    archivoSalida << "\n* Matriz:" << endl;
    for (const auto &fila : resultado.matrizScores) {
      for (size_t j = 0; j < fila.size(); ++j) {
        archivoSalida << setw(4) << fila[j] << (j == fila.size() - 1 ? "" : "\t");
      }
      archivoSalida << endl;
    }

    archivoSalida << "\n* N de alineamientos optimos: " << resultado.alineamientos.size() << endl;

    archivoSalida << "\n* Alineamientos optimos:" << endl;
    for (size_t k = 0; k < resultado.alineamientos.size(); ++k) {
      const auto &info = resultado.alineamientos[k];
      archivoSalida << "\t*Alineamiento " << k + 1 << ":" << endl;
      // subsecuencia comun
      archivoSalida << "\t\tS1: " << info.s1_alineada << endl;
      archivoSalida << "\t\tS2: " << info.s2_alineada << endl;
      // posición donde se encuentran ambas cadenas
      archivoSalida << "\t\tS1_pos: [" << info.start_s1 << " - " << info.end_s1 << "]" << endl;
      archivoSalida << "\t\tS2_pos: [" << info.start_s2 << " - " << info.end_s2 << "]" << endl << endl;
    }

    archivoSalida.close();
    cout << "Resultados guardados en " << nombreArchivo << endl;
  } else {
    cerr << "Error al abrir el archivo " << nombreArchivo << endl;
  }
}

int main() {
  cout << "--- Alineamiento Local (Smith-Waterman) ---" << endl;

  string sec1 = "AGCT";
  string sec2 = "GCA";
  cout << "\nProcesando S1: " << sec1 << " y S2: " << sec2 << endl;
  ResultadoAlineamientoLocal res1 = alineamientoLocal(sec1, sec2);
  guardarResultados("resultado_alineamiento_1.txt", res1);

  string sec3 = "CCCGGGTTTAAA";
  string sec4 = "TTTGGGCCCAAA";
  cout << "\nProcesando S3: " << sec3 << " y S4: " << sec4 << endl;
  ResultadoAlineamientoLocal res2 = alineamientoLocal(sec3, sec4);
  guardarResultados("resultado_alineamiento_2.txt", res2);

  string sec5 = "GGTTGACTA";
  string sec6 = "TGTTAGGG";
  cout << "\nProcesando S5: " << sec5 << " y S6: " << sec6 << endl;
  ResultadoAlineamientoLocal res3 = alineamientoLocal(sec5, sec6);
  guardarResultados("resultado_alineamiento_3.txt", res3);

  return 0;
}
