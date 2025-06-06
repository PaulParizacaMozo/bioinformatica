#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// Params alineamiento global
const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

// Estructura para el resultado de un alineamiento par-a-par global
struct ResultadoAlineamientoPar {
  string sec1Alineada;
  string sec2Alineada;
  int score;
};

// Implementacion alineamiento global
ResultadoAlineamientoPar alineamientoGlobalPar(const string &sec1, const string &sec2) {
  int longitud1 = sec1.length();
  int longitud2 = sec2.length();
  vector<vector<int>> matriz(longitud1 + 1, vector<int>(longitud2 + 1));

  for (int i = 0; i <= longitud1; ++i)
    matriz[i][0] = i * GAP;
  for (int j = 0; j <= longitud2; ++j)
    matriz[0][j] = j * GAP;

  for (int i = 1; i <= longitud1; ++i) {
    for (int j = 1; j <= longitud2; ++j) {
      int sumaResta = (sec1[i - 1] == sec2[j - 1] ? MATCH : MISMATCH);
      int scoreDiagonal = matriz[i - 1][j - 1] + sumaResta;
      int scoreArriba = matriz[i - 1][j] + GAP;
      int scoreIzquierda = matriz[i][j - 1] + GAP;
      matriz[i][j] = max({scoreDiagonal, scoreArriba, scoreIzquierda});
    }
  }

  ResultadoAlineamientoPar res;
  res.score = matriz[longitud1][longitud2];
  string alineado1_rev = "", alineado2_rev = "";
  int i = longitud1, j = longitud2;

  while (i > 0 || j > 0) {
    // Caso Diagonal
    if (i > 0 && j > 0 && matriz[i][j] == matriz[i - 1][j - 1] + (sec1[i - 1] == sec2[j - 1] ? MATCH : MISMATCH)) {
      alineado1_rev += sec1[i - 1];
      alineado2_rev += sec2[j - 1];
      i--;
      j--;
    }
    // Caso arriba
    else if (i > 0 && matriz[i][j] == matriz[i - 1][j] + GAP) {
      alineado1_rev += sec1[i - 1];
      alineado2_rev += '-';
      i--;
    }
    // Caso izquierda
    else if (j > 0 && matriz[i][j] == matriz[i][j - 1] + GAP) {
      alineado1_rev += '-';
      alineado2_rev += sec2[j - 1];
      j--;
    }
    // Caso limite cuando una sec se acaba antes
    else {
      if (i > 0) {
        alineado1_rev += sec1[i - 1];
        alineado2_rev += '-';
        i--;
      } else if (j > 0) {
        alineado1_rev += '-';
        alineado2_rev += sec2[j - 1];
        j--;
      } else { // Ambos i y j son 0
        break;
      }
    }
  }
  reverse(alineado1_rev.begin(), alineado1_rev.end());
  reverse(alineado2_rev.begin(), alineado2_rev.end());
  res.sec1Alineada = alineado1_rev;
  res.sec2Alineada = alineado2_rev;
  return res;
}

// Estructura para el resultado del Alineamiento Estrella
struct ResultadoAlineamientoEstrella {
  vector<vector<int>> matrizScores;
  int indiceSecCentralOriginal;
  vector<ResultadoAlineamientoPar> alineamientosConEstrella;
  vector<string> alineamientoMultiple;
};

// Implementación del Alineamiento Estrella
ResultadoAlineamientoEstrella alineamientoEstrella(const vector<string> &secs) {
  ResultadoAlineamientoEstrella resultado;
  int numsecs = secs.size();

  if (numsecs < 2) {
    if (numsecs == 1) {
      resultado.indiceSecCentralOriginal = 0;
      resultado.alineamientoMultiple.push_back(secs[0]);
    }
    return resultado;
  }

  // 1. Calcular la matriz de scores
  resultado.matrizScores.assign(numsecs, vector<int>(numsecs, 0));
  vector<long long> sumaScores(numsecs, 0);

  for (int i = 0; i < numsecs; ++i) {
    for (int j = i + 1; j < numsecs; ++j) {
      // Realizamos un alineamiento global entre las secuencias i y j
      ResultadoAlineamientoPar resAlineamiento = alineamientoGlobalPar(secs[i], secs[j]);
      // Guardamos el score en la matriz de scores, que es simétrica
      resultado.matrizScores[i][j] = resultado.matrizScores[j][i] = resAlineamiento.score;
      // Sumamos los scores de cada secuencia para determinar cuál tendrá la mayor relación con las otras
      sumaScores[i] += resAlineamiento.score;
      sumaScores[j] += resAlineamiento.score;
    }
  }

  // 2. Seleccionar la secuencia estrella
  resultado.indiceSecCentralOriginal = 0; // Inicializamos el índice de la secuencia central
  long long scoreTotalMaximo = -3e18;     // Usamos un valor muy pequeño para iniciar la comparación

  for (int i = 0; i < numsecs; ++i) {
    // Buscamos la secuencia con la mayor suma de scores
    if (sumaScores[i] > scoreTotalMaximo) {
      scoreTotalMaximo = sumaScores[i];       // Actualizamos el score máximo
      resultado.indiceSecCentralOriginal = i; // Guardamos el índice de la secuencia estrella
    }
  }

  // Guardamos la secuencia original que corresponde a la estrella
  string secCentralOriginalStr = secs[resultado.indiceSecCentralOriginal];

  // 3. Alinear todas las otras secuencias con la estrella
  map<int, int> mapaIndiceOriginalAIndiceAlineamientoParApar; // Mapea el índice original al índice en los alineamientos
  for (int i = 0; i < numsecs; ++i) {
    if (i != resultado.indiceSecCentralOriginal) {
      // Guardamos los alineamientos de las secuencias con la secuencia estrella
      mapaIndiceOriginalAIndiceAlineamientoParApar[i] = resultado.alineamientosConEstrella.size();
      // Realizamos un alineamiento global entre la secuencia estrella y la secuencia i
      resultado.alineamientosConEstrella.push_back(alineamientoGlobalPar(secCentralOriginalStr, secs[i]));
    }
  }

  // 4. Construir Alineamiento Múltiple (MSA)
  vector<string> filasMSAFinal(numsecs);
  vector<int> punterosAlineamientoParApar(numsecs, 0);
  int punteroSecCentralOriginal = 0;

  while (true) {
    // Primero, verificamos si hemos terminado
    bool secCentralTerminado = (punteroSecCentralOriginal >= secCentralOriginalStr.length());
    bool todosLosAlineamientosTerminados = true;
    for (auto const &[indiceOriginal, indiceAlineamiento] : mapaIndiceOriginalAIndiceAlineamientoParApar) {
      if (punterosAlineamientoParApar[indiceOriginal] <
          resultado.alineamientosConEstrella[indiceAlineamiento].sec1Alineada.length()) {
        todosLosAlineamientosTerminados = false;
        break;
      }
    }
    if (secCentralTerminado && todosLosAlineamientosTerminados)
      break;

    // Paso 1: Buscar si alguna secuencia tiene una inserción en la posición actual
    bool insercionProcesada = false;
    for (auto const &[indiceOriginal, indiceAlineamiento] : mapaIndiceOriginalAIndiceAlineamientoParApar) {
      const auto &alineamiento = resultado.alineamientosConEstrella[indiceAlineamiento];
      int &punteroActual = punterosAlineamientoParApar[indiceOriginal];

      if (punteroActual < alineamiento.sec1Alineada.length() && alineamiento.sec1Alineada[punteroActual] == '-') {
        insercionProcesada = true;
        // La secuencia 'indiceOriginal' tiene una inserción. Crear una columna en el MSA.
        for (int i = 0; i < numsecs; ++i) {
          if (i == indiceOriginal) {
            filasMSAFinal[i] += alineamiento.sec2Alineada[punteroActual];
          } else {
            filasMSAFinal[i] += '-';
          }
        }
        punteroActual++;
      }
    }

    if (insercionProcesada)
      continue; // Volver a empezar para procesar todas las inserciones en este punto

    // Paso 2: Si no hubo inserciones, procesar la columna de alineamiento normal
    filasMSAFinal[resultado.indiceSecCentralOriginal] += secCentralOriginalStr[punteroSecCentralOriginal];

    for (auto const &[indiceOriginal, indiceAlineamiento] : mapaIndiceOriginalAIndiceAlineamientoParApar) {
      const auto &alineamiento = resultado.alineamientosConEstrella[indiceAlineamiento];
      int &punteroActual = punterosAlineamientoParApar[indiceOriginal];

      if (punteroActual < alineamiento.sec2Alineada.length()) {
        filasMSAFinal[indiceOriginal] += alineamiento.sec2Alineada[punteroActual];
      } else {
        filasMSAFinal[indiceOriginal] += '-';
      }
      punteroActual++;
    }
    punteroSecCentralOriginal++;
  }
  resultado.alineamientoMultiple = filasMSAFinal;

  return resultado;
}

// Función para guardar los resultados del Alineamiento Estrella
void guardarResultadosAlineamientoEstrella(const string &nombreArchivo, const ResultadoAlineamientoEstrella &resultado,
                                           const vector<string> &secs) {
  ofstream archivoSalida(nombreArchivo);
  if (!archivoSalida.is_open()) {
    cerr << "Error al abrir el archivo " << nombreArchivo << endl;
    return;
  }

  archivoSalida << "* Matriz de scores (alineamientos globales par-a-par):" << endl;
  if (!resultado.matrizScores.empty() && !resultado.matrizScores[0].empty()) {
    for (size_t i = 0; i < resultado.matrizScores.size(); ++i) {
      for (size_t j = 0; j < resultado.matrizScores[i].size(); ++j) {
        archivoSalida << setw(5) << resultado.matrizScores[i][j];
      }
      archivoSalida << endl;
    }
  } else {
    archivoSalida << "No disponible (por ejemplo, si solo hay 1 sec o ninguna)." << endl;
  }

  if (secs.empty()) {
    archivoSalida << "\nNo hay secs para procesar." << endl;
  } else {
    archivoSalida << "\nsec estrella (índice original): " << resultado.indiceSecCentralOriginal << " ("
                  << secs[resultado.indiceSecCentralOriginal] << ")" << endl;

    archivoSalida << "\n* Alineamientos de cada sec con la estrella:" << endl;
    int conteosecsProcesadas = 0;
    for (size_t i = 0; i < secs.size(); ++i) {
      if ((int)i == resultado.indiceSecCentralOriginal)
        continue; // Saltar la propia estrella

      // Encontrar el alineamiento correspondiente en resultado.alineamientosConEstrella
      if (conteosecsProcesadas < resultado.alineamientosConEstrella.size()) {
        const auto &parejaAlineada = resultado.alineamientosConEstrella[conteosecsProcesadas];
        archivoSalida << "Alineamiento con S" << resultado.indiceSecCentralOriginal << " (estrella) y S" << i << " (" << secs[i]
                      << "):" << endl;
        archivoSalida << "  Estrella: " << parejaAlineada.sec1Alineada << endl;
        archivoSalida << "  Sec " << i << "    : " << parejaAlineada.sec2Alineada << endl;
        archivoSalida << "  score : " << parejaAlineada.score << endl << endl;
        conteosecsProcesadas++;
      }
    }

    archivoSalida << "\n* Alineamiento múltiple:" << endl;
    if (!resultado.alineamientoMultiple.empty()) {
      for (size_t i = 0; i < resultado.alineamientoMultiple.size(); ++i) {
        archivoSalida << "S" << i << " (" << secs[i] << "): " << resultado.alineamientoMultiple[i] << endl;
      }
    } else {
      archivoSalida << "No disponible." << endl;
    }
  }

  archivoSalida.close();
  cout << "Resultados de Alineamiento Estrella guardados en " << nombreArchivo << endl;
}

int main() {
  cout << "--- Laboratorio de Alineamiento Estrella ---" << endl;
  vector<std::string> secs1 = {"ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT", "ACTGACC"};

  cout << "\nProcesando Conjunto 1 de secs..." << endl;
  ResultadoAlineamientoEstrella resAlineamiento1 = alineamientoEstrella(secs1);
  guardarResultadosAlineamientoEstrella("alineamiento_estrella_1.txt", resAlineamiento1, secs1);

  vector<string> secs2 = {"MWAFGGRAAVGLLPRTASRASAWVGNPRWREPIVTCGRRGLHVTVNAGATRHAHLNLHYLQILNIKKQSV",
                          "CVVHLRNLGTLDNPSSLDETAYERLAEETLDSLAEFFEDLADKPYTLEDYDVSFGDGVLTIKLGGDLGTY",
                          "VINKQTPNKQIWLSSPSSGPKRYDWTGKNWVYSHDGVSLHELLARELTKALNTKLDLSSLAYSGKGT"};

  cout << "\nProcesando Conjunto 2 de secs..." << endl;
  ResultadoAlineamientoEstrella resAlineamiento2 = alineamientoEstrella(secs2);
  guardarResultadosAlineamientoEstrella("alineamiento_estrella_2.txt", resAlineamiento2, secs2);
  return 0;
}
