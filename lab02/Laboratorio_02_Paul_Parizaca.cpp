#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

// Función para verificar si una cadena es subcadena de otra
bool esSubstring(const string &cadena, const string &subcadena) {
  if (subcadena.length() > cadena.length()) {
    return false;
  }
  for (size_t i = 0; i <= cadena.length() - subcadena.length(); ++i) {
    bool coincide = true;
    for (size_t j = 0; j < subcadena.length(); ++j) {
      if (cadena[i + j] != subcadena[j]) {
        coincide = false;
        break;
      }
    }
    if (coincide) {
      return true;
    }
  }
  return false;
}

// Función para calcular el score entre dos cadenas
int calcularScoreSimple(const string &cadena1, const string &cadena2) {
  int score = 0;
  int longitudMinima = min(cadena1.length(), cadena2.length());
  for (int i = 0; i < longitudMinima; ++i) {
    if (cadena1[i] == cadena2[i]) {
      score += 1;
    } else {
      score -= 2;
    }
  }
  // Penalizar por la diferencia de longitud
  score -= 2 * abs((int)cadena1.length() - (int)cadena2.length());
  return score;
}

// Parámetros para el alineamiento global
const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

// Estructura para los resultados del alineamiento
struct ResultadoAlineamiento {
  int scoreFinal;
  vector<vector<int>> matrizScores;
  int cantidadAlineamientos;
  vector<pair<string, string>> alineamientosGenerados;
};

// Imprimir matriz
void imprimirMatriz(const vector<vector<int>> &matriz) {
  for (const auto &fila : matriz) {
    for (int val : fila) {
      cout << setw(4) << val << " ";
    }
    cout << endl;
  }
}

// Función reconstruir para encontrar todos los alineamientos óptimos
void reconstruir(const string &s1, const string &s2, const vector<vector<int>> &matriz, int i, int j, string alin1,
                 string alin2, vector<pair<string, string>> &alineamientos) {
  if (i == 0 && j == 0) {
    reverse(alin1.begin(), alin1.end());
    reverse(alin2.begin(), alin2.end());
    alineamientos.push_back({alin1, alin2});
    return;
  }

  int scoreActual = matriz[i][j];
  int scoreDiagonal = (i > 0 && j > 0) ? matriz[i - 1][j - 1] : -1e9; // Un valor muy pequeño si está fuera de límites
  int scoreArriba = (i > 0) ? matriz[i - 1][j] : -1e9;
  int scoreIzquierda = (j > 0) ? matriz[i][j - 1] : -1e9;

  // Camino diagonal
  if (i > 0 && j > 0 && scoreActual == scoreDiagonal + (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH)) {
    reconstruir(s1, s2, matriz, i - 1, j - 1, alin1 + s1[i - 1], alin2 + s2[j - 1], alineamientos);
  }

  // Camino desde arriba
  if (i > 0 && scoreActual == scoreArriba + GAP) {
    reconstruir(s1, s2, matriz, i - 1, j, alin1 + s1[i - 1], alin2 + '-', alineamientos);
  }

  // Camino desde la izquierda
  if (j > 0 && scoreActual == scoreIzquierda + GAP) {
    reconstruir(s1, s2, matriz, i, j - 1, alin1 + '-', alin2 + s2[j - 1], alineamientos);
  }
}

// Implementación del alineamiento global (Needleman-Wunch)
ResultadoAlineamiento alineamientoGlobal(const string &s1, const string &s2) {
  int n = s1.length();
  int m = s2.length();

  vector<vector<int>> matriz(n + 1, vector<int>(m + 1));

  // Inicializar la matriz de scores
  for (int i = 0; i <= n; ++i) {
    matriz[i][0] = i * GAP;
  }
  for (int j = 0; j <= m; ++j) {
    matriz[0][j] = j * GAP;
  }

  // Llenar la matriz de scores
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int scoreDiagonal = matriz[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH);
      int scoreArriba = matriz[i - 1][j] + GAP;
      int scoreIzquierda = matriz[i][j - 1] + GAP;
      matriz[i][j] = max({scoreDiagonal, scoreArriba, scoreIzquierda});
    }
  }

  ResultadoAlineamiento resultado;
  resultado.scoreFinal = matriz[n][m];
  resultado.matrizScores = matriz;

  // Realizar reconstruccion para obtener los alineamientos
  vector<pair<string, string>> alineamientos;
  reconstruir(s1, s2, matriz, n, m, "", "", alineamientos);
  resultado.alineamientosGenerados = alineamientos;
  resultado.cantidadAlineamientos = alineamientos.size();

  return resultado;
}

// Función para guardar resultados
void guardarResultados(const string &nombreArchivo, const ResultadoAlineamiento &resultado) {
  ofstream archivoSalida(nombreArchivo);
  if (archivoSalida.is_open()) {
    archivoSalida << "* Score final(Optimo): " << resultado.scoreFinal << endl;
    archivoSalida << "\n* Matriz:" << endl;
    for (const auto &fila : resultado.matrizScores) {
      for (size_t j = 0; j < fila.size(); ++j) {
        archivoSalida << setw(4) << fila[j] << (j == fila.size() - 1 ? "" : "\t");
      }
      archivoSalida << endl;
    }

    archivoSalida << "\n* N de alineamientos optimos: " << resultado.cantidadAlineamientos << endl;
    archivoSalida << "\n* Alineamientos optimos:" << endl;
    for (size_t i = 0; i < resultado.alineamientosGenerados.size(); ++i) {
      archivoSalida << "\t*Alineamiento " << i + 1 << ":" << endl;
      archivoSalida << "\t\t" << resultado.alineamientosGenerados[i].first << endl;
      archivoSalida << "\t\t" << resultado.alineamientosGenerados[i].second << endl << endl;
    }

    archivoSalida.close();
    cout << "Resultados guardados en " << nombreArchivo << endl;
  } else {
    cerr << "Error al abrir el archivo " << nombreArchivo << endl;
  }
}

int main() {
  // Cadenas con diferentes tamaños (larga, mediana y corta)
  string secA = "GCTAGGCGATCGGCTAAGGCTAGTACGATGCA";
  string secB = "CGATCGGCTAAGGCTAGT";
  string secC = "TACG";
  string secD = "GCTAGGCGATCGGCTAAGGC";
  string secE = "CGATCGGCTAAGGC";

  // 1. Substring
  cout << "--- Substring ---" << endl;
  cout << "'" << secC << "' es substring de '" << secA << "': " << (esSubstring(secA, secC) ? "Si" : "No") << endl;
  cout << "'" << secB << "' es substring de '" << secA << "': " << (esSubstring(secA, secB) ? "Si" : "No") << endl;
  cout << "'" << secC << "' es substring de '" << secB << "': " << (esSubstring(secB, secC) ? "Si" : "No") << endl;
  cout << "'" << secC << "' es substring de '" << secD << "': " << (esSubstring(secD, secA) ? "Si" : "No") << endl;
  cout << endl;

  // 2. Calcular Score Simple
  cout << "--- Calculo de Score Simple ---" << endl;
  cout << "Score entre '" << secA << "' y '" << secA << "': " << calcularScoreSimple(secA, secA) << endl;
  cout << "Score entre '" << secE << "' y '" << secB << "': " << calcularScoreSimple(secE, secB) << endl;
  cout << "Score entre '" << secA << "' y '" << secD << "': " << calcularScoreSimple(secA, secD) << endl;
  cout << "Score entre '" << secB << "' y '" << secC << "': " << calcularScoreSimple(secB, secC) << endl;
  cout << endl;

  // 3. Alineamiento Global
  cout << "--- Alineamiento Global ---" << endl;
  string sec1 = "GATTACA";
  string sec2 = "GCATGCU";
  cout << "\nAlineando '" << sec1 << "' y '" << sec2 << "'" << endl;
  ResultadoAlineamiento resAG = alineamientoGlobal(sec1, sec2);
  guardarResultados("alineamiento_global_1.txt", resAG);

  string sec3 = "ATGCGTACG";
  string sec4 = "GCTAGC";
  cout << "\nAlineando '" << sec3 << "' y '" << sec4 << "'" << endl;
  ResultadoAlineamiento resAG2 = alineamientoGlobal(sec3, sec4);
  guardarResultados("alineamiento_global_2.txt", resAG2);

  string sec5 = "CGTAGCTAGCTACGAT";
  string sec6 = "AGCTGACTG";
  cout << "\nAlineando '" << sec5 << "' y '" << sec6 << "'" << endl;
  ResultadoAlineamiento resAG3 = alineamientoGlobal(sec5, sec6);
  guardarResultados("alineamiento_global_3.txt", resAG3);

  return 0;
}
