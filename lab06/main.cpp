#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Cluster {
  int id;
  string etiqueta;
  int tam;
  bool activo;
};

using MatrizEtiquetas = pair<vector<vector<double>>, vector<string>>;

MatrizEtiquetas leerMatriz(const string &rutaArchivo, bool conEtiquetas) {
  ifstream archivo(rutaArchivo);
  if (!archivo.is_open()) {
    cerr << "Error: No se pudo abrir el archivo " << rutaArchivo << endl;
    exit(1);
  }
  string linea;
  vector<string> etiquetas;
  int n;

  if (conEtiquetas) {
    getline(archivo, linea);
    stringstream ss(linea);
    string etiqueta;
    while (ss >> etiqueta)
      etiquetas.push_back(etiqueta);
    n = etiquetas.size();
  }

  vector<vector<double>> datos;
  while (getline(archivo, linea)) {
    if (linea.empty() || linea.find_first_not_of(" \t\n\v\f\r") == string::npos)
      continue;
    datos.emplace_back();
    stringstream ss(linea);
    double valor;
    while (ss >> valor)
      datos.back().push_back(valor);
  }

  if (!conEtiquetas) {
    n = datos.size() + 1;
    for (int i = 0; i < n; ++i)
      etiquetas.push_back("S" + to_string(i + 1));
  }

  if (datos.size() != n - 1) {
    cerr << "Error: Conflicto de dimensiones en la matriz." << endl;
    exit(1);
  }

  vector<vector<double>> matriz(n, vector<double>(n, 0.0));
  for (int i = 0; i < n - 1; ++i) {
    for (int j = 0; j < datos[i].size(); ++j) {
      matriz[i + 1][j] = datos[i][j];
      matriz[j][i + 1] = datos[i][j];
    }
  }
  return {matriz, etiquetas};
}

string formatearMatriz(const vector<vector<double>> &matriz, const vector<Cluster> &clusters) {
  stringstream ss;
  vector<string> etiquetas_activas;
  vector<int> indices_activos;

  for (int i = 0; i < clusters.size(); ++i) {
    if (clusters[i].activo) {
      etiquetas_activas.push_back(clusters[i].etiqueta);
      indices_activos.push_back(i);
    }
  }

  ss << "\t";
  for (const auto &etiqueta : etiquetas_activas) {
    ss << etiqueta << "\t";
  }
  ss << "\n";

  for (size_t i = 0; i < indices_activos.size(); ++i) {
    ss << etiquetas_activas[i] << "\t";
    for (size_t j = 0; j < indices_activos.size(); ++j) {
      ss << fixed << setprecision(2) << matriz[indices_activos[i]][indices_activos[j]] << "\t";
    }
    ss << "\n";
  }
  return ss.str();
}

void clusteringYGenerarSalidas(const string &archivoLog, const string &archivoEnlace, vector<vector<double>> &matriz,
                               const vector<string> &etiquetasIniciales, const string &metodoEnlace) {

  ofstream logStream(archivoLog);
  ofstream enlaceStream(archivoEnlace);

  logStream << "--- INICIO CLUSTERING CON METODO: " << metodoEnlace.c_str() << " ---\n\n";

  int n = etiquetasIniciales.size();
  vector<Cluster> clusters;
  for (int i = 0; i < n; ++i) {
    clusters.push_back({i, etiquetasIniciales[i], 1, true});
  }

  int proximoId = n;
  int clustersActivos = n;
  int paso = 1;

  while (clustersActivos > 1) {
    logStream << "PASO " << paso++ << ":\n";
    logStream << "Matriz de Distancias Actual:\n";
    logStream << formatearMatriz(matriz, clusters) << "\n";

    double distMin = numeric_limits<double>::max();
    int idx1 = -1, idx2 = -1;

    for (int i = 0; i < n; ++i) {
      if (!clusters[i].activo)
        continue;
      for (int j = i + 1; j < n; ++j) {
        if (!clusters[j].activo)
          continue;
        if (matriz[i][j] < distMin) {
          distMin = matriz[i][j];
          idx1 = i;
          idx2 = j;
        }
      }
    }

    logStream << "-> Se unen los clusters '" << clusters[idx1].etiqueta << "' y '" << clusters[idx2].etiqueta << "'.\n";
    logStream << "-> Distancia de union: " << fixed << setprecision(2) << distMin << "\n\n";

    int nuevoTam = clusters[idx1].tam + clusters[idx2].tam;
    enlaceStream << fixed << setprecision(4) << clusters[idx1].id << " " << clusters[idx2].id << " " << distMin << " "
                 << nuevoTam << "\n";

    for (int k = 0; k < n; ++k) {
      if (!clusters[k].activo || k == idx1 || k == idx2)
        continue;
      double nuevaDist;
      if (metodoEnlace == "minima")
        nuevaDist = min(matriz[idx1][k], matriz[idx2][k]);
      else if (metodoEnlace == "maxima")
        nuevaDist = max(matriz[idx1][k], matriz[idx2][k]);
      else
        nuevaDist = (matriz[idx1][k] * clusters[idx1].tam + matriz[idx2][k] * clusters[idx2].tam) / nuevoTam;
      matriz[idx1][k] = matriz[k][idx1] = nuevaDist;
    }

    clusters[idx1].id = proximoId;
    clusters[idx1].etiqueta = "(" + clusters[idx1].etiqueta + "," + clusters[idx2].etiqueta + ")";
    clusters[idx1].tam = nuevoTam;

    clusters[idx2].activo = false;
    clustersActivos--;
    proximoId++;
  }

  logStream << "--- CLUSTERING FINALIZADO ---\n\n";
  logStream << "Dendrograma (formato Newick):\n";
  for (const auto &c : clusters) {
    if (c.activo) {
      logStream << c.etiqueta << "\n";
      break;
    }
  }

  logStream.close();
  enlaceStream.close();
}

int main() {
  cout << "--- Cluster Calculator ---\n";

  string rutaArchivo;
  cout << "Introduce la ruta al archivo .txt con la matriz: ";
  cin >> rutaArchivo;

  int opcion;
  cout << "\nFormato de la matriz:\n1. Triangular (sin etiquetas)\n2. Triangular (con etiquetas)\nElige: ";
  cin >> opcion;

  MatrizEtiquetas datos = (opcion == 2) ? leerMatriz(rutaArchivo, true) : leerMatriz(rutaArchivo, false);

  string etiquetasStr;
  for (size_t i = 0; i < datos.second.size(); ++i) {
    etiquetasStr += datos.second[i] + (i == datos.second.size() - 1 ? "" : ",");
  }

  vector<string> metodos = {"minima", "maxima", "promedio"};

  for (const auto &metodo : metodos) {
    cout << "\nProcesando metodo: '" << metodo << "'...\n";

    vector<vector<double>> matrizCopia = datos.first;
    string archivoLog = "log_" + metodo + ".txt";
    string archivoEnlace = "enlace_" + metodo + ".txt";

    clusteringYGenerarSalidas(archivoLog, archivoEnlace, matrizCopia, datos.second, metodo);
    cout << "-> Archivo de log detallado generado: '" << archivoLog << "'\n";
    cout << "-> Archivo de enlace para Python generado: '" << archivoEnlace << "'\n";

    string archivoPng = "dendrograma_" + metodo + ".png";
    string comando =
        "python plotter.py \"" + archivoEnlace + "\" \"" + etiquetasStr + "\" \"" + archivoPng + "\" \"" + metodo + "\"";

    cout << "-> Ejecutando comando: " << comando << "\n";
    int resultado = system(comando.c_str());

    if (resultado == 0) {
      cout << "-> Grafico guardado en '" << archivoPng << "'\n";
    } else {
      cerr << "-> Error al ejecutar el script de Python.\n";
    }
  }

  cout << "\nProceso terminado\n";
  return 0;
}
