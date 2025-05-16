#include <algorithm> // Para transform (convertir a mayúsculas)
#include <fstream>   // Para manejo de archivos (ifstream)
#include <iostream>  // Para entrada y salida estándar (cout, cerr)
#include <map>       // Para usar map (tablas de codones y nombres de aminoácidos)
#include <string>    // Para usar la clase string
using namespace std;

// Constante global: Tabla de traducción de codones de ARN a nombres de aminoácidos.
// Un codón es una secuencia de tres nucleótidos.
const map<string, string> tablaCodonesARN = {
    {"UUU", "Fenilalanina"},    {"UUC", "Fenilalanina"},    {"UUA", "Leucina"},         {"UUG", "Leucina"},         {"CUU", "Leucina"},
    {"CUC", "Leucina"},         {"CUA", "Leucina"},         {"CUG", "Leucina"},         {"AUU", "Isoleucina"},      {"AUC", "Isoleucina"},
    {"AUA", "Isoleucina"},      {"AUG", "Metionina"},       {"GUU", "Valina"},          {"GUC", "Valina"},          {"GUA", "Valina"},
    {"GUG", "Valina"},          {"UCU", "Serina"},          {"UCC", "Serina"},          {"UCA", "Serina"},          {"UCG", "Serina"},
    {"CCU", "Prolina"},         {"CCC", "Prolina"},         {"CCA", "Prolina"},         {"CCG", "Prolina"},         {"ACU", "Treonina"},
    {"ACC", "Treonina"},        {"ACA", "Treonina"},        {"ACG", "Treonina"},        {"GCU", "Alanina"},         {"GCC", "Alanina"},
    {"GCA", "Alanina"},         {"GCG", "Alanina"},         {"UAU", "Tirosina"},        {"UAC", "Tirosina"},        {"UAA", "DETENCION"},
    {"UAG", "DETENCION"},       {"UGA", "DETENCION"},       {"CAU", "Histidina"},       {"CAC", "Histidina"},       {"CAA", "Glutamina"},
    {"CAG", "Glutamina"},       {"AAU", "Asparagina"},      {"AAC", "Asparagina"},      {"AAA", "Lisina"},          {"AAG", "Lisina"},
    {"GAU", "Acido Aspartico"}, {"GAC", "Acido Aspartico"}, {"GAA", "Acido Glutamico"}, {"GAG", "Acido Glutamico"}, {"UGU", "Cisteina"},
    {"UGC", "Cisteina"},        {"UGG", "Triptofano"},      {"CGU", "Arginina"},        {"CGC", "Arginina"},        {"CGA", "Arginina"},
    {"CGG", "Arginina"},        {"AGU", "Serina"},          {"AGC", "Serina"},          {"AGA", "Arginina"},        {"AGG", "Arginina"},
    {"GGU", "Glicina"},         {"GGC", "Glicina"},         {"GGA", "Glicina"},         {"GGG", "Glicina"}};

// Constante global: Mapa para convertir el código de una letra de un aminoácido a su nombre completo.
// SLA = Single Letter Aminoacid (Aminoácido de una Sola Letra)
const map<char, string> nombresAminoacidosSLA = {{'A', "Alanina"},   {'R', "Arginina"},     {'N', "Asparagina"},      {'D', "Acido Aspartico"},
                                                 {'C', "Cisteina"},  {'Q', "Glutamina"},    {'E', "Acido Glutamico"}, {'G', "Glicina"},
                                                 {'H', "Histidina"}, {'I', "Isoleucina"},   {'L', "Leucina"},         {'K', "Lisina"},
                                                 {'M', "Metionina"}, {'F', "Fenilalanina"}, {'P', "Prolina"},         {'S', "Serina"},
                                                 {'T', "Treonina"},  {'V', "Valina"},       {'W', "Triptofano"},      {'Y', "Tirosina"}};

// Función: obtenerAminoacidoDeCodon
// Propósito: Traduce un codón (de ADN o ARN) al aminoácido correspondiente o señal de detención.
// Parámetros:
//   - codon: La secuencia de 3 bases a traducir.
//   - esADN: Booleano que indica si el codón es de ADN (true) o ARN (false). Si es ADN, 'T' se convierte a 'U'.
// Retorna: El nombre del aminoácido o "DETENCION", o un mensaje de error si el codón es inválido.
string obtenerAminoacidoDeCodon(const string &codon, bool esADN) {
  if (codon.length() != 3) { // Un codón siempre tiene 3 bases
    return "Longitud de codon invalida";
  }
  string codonARN = codon;
  // if (esADN) {
  //   // Si es ADN, convertir Timina (T) a Uracilo (U) para usar la tabla de codones de ARN
  //   transform(codonARN.begin(), codonARN.end(), codonARN.begin(), [](unsigned char c) { return (c == 'T') ? 'U' : c; });
  // }

  // Buscar el codón de ARN en la tabla
  auto iterador = tablaCodonesARN.find(codonARN);
  if (iterador != tablaCodonesARN.end()) {
    return iterador->second; // Retorna el nombre del aminoácido o "DETENCION"
  }
  return "Codon desconocido"; // Si el codón no se encuentra en la tabla
}

// Función: obtenerNombreCompletoAminoacido
// Propósito: Obtiene el nombre completo de un aminoácido a partir de su código de una letra (SLA).
// Parámetros:
//   - sla: El carácter que representa el código de una letra del aminoácido.
// Retorna: El nombre completo del aminoácido o "Aminoacido desconocido".
string obtenerNombreCompletoAminoacido(char sla) {
  auto iterador = nombresAminoacidosSLA.find(sla);
  if (iterador != nombresAminoacidosSLA.end()) {
    return iterador->second; // Retorna el nombre completo
  }
  return "Aminoacido desconocido"; // Si el código SLA no se encuentra
}

// Función: procesarSecuencia
// Propósito: Analiza una secuencia de entrada para determinar si es ADN, ARN o proteína,
//            e imprime la clasificación y detalles adicionales.
// Parámetros:
//   - secuencia_original: La cadena de texto leída del archivo.
void procesarSecuencia(const string &secuencia_original) {
  string secuencia_procesada = secuencia_original;
  // Convertir toda la secuencia a mayúsculas para un análisis no sensible a mayúsculas/minúsculas
  transform(secuencia_procesada.begin(), secuencia_procesada.end(), secuencia_procesada.begin(), ::toupper);

  // cout << "Secuencia original: \"" << secuencia_original << "\" -> ";

  if (secuencia_procesada.empty()) { // Comprobar si la secuencia está vacía
    cout << "Vacia." << endl;
    return; // No hay más que procesar para una secuencia vacía
  }

  // Banderas (flags) para caracterizar la secuencia
  bool contiene_T = false;                  // Indica si la secuencia contiene Timina ('T')
  bool contiene_U = false;                  // Indica si la secuencia contiene Uracilo ('U')
  bool solo_caracteres_ACGT = true;         // Asume inicialmente que solo contiene bases de ADN (A,C,G,T)
  bool solo_caracteres_ACGU = true;         // Asume inicialmente que solo contiene bases de ARN (A,C,G,U)
  bool solo_caracteres_SLA_proteina = true; // Asume inicialmente que solo contiene códigos de aminoácidos de una letra
  bool contiene_caracter_invalido = false;  // Indica si se encontró algún carácter no permitido

  // Conjuntos de caracteres válidos para cada tipo de molécula
  const string conjunto_SLA_proteina = "ACDEFGHIKLMNPQRSTVWY"; // Códigos válidos de una letra para aminoácidos
  const string bases_ADN = "ACGT";                             // Bases canónicas del ADN
  const string bases_ARN = "ACGU";                             // Bases canónicas del ARN

  // Iterar sobre cada carácter de la secuencia para evaluarlo
  for (char caracter_actual : secuencia_procesada) {
    if (caracter_actual == 'T') {
      contiene_T = true;
    } else if (caracter_actual == 'U') {
      contiene_U = true;
    }

    if (bases_ADN.find(caracter_actual) == string::npos) {
      solo_caracteres_ACGT = false;
    }
    if (bases_ARN.find(caracter_actual) == string::npos) {
      solo_caracteres_ACGU = false;
    }
    if (conjunto_SLA_proteina.find(caracter_actual) == string::npos) {
      solo_caracteres_SLA_proteina = false;
    }

    if (caracter_actual != 'U' && conjunto_SLA_proteina.find(caracter_actual) == string::npos) {
      contiene_caracter_invalido = true;
      break;
    }
  }

  // Lógica de decisión basada en las banderas evaluadas
  if (contiene_caracter_invalido) {
    cout << "Contiene caracteres no validos." << endl;
  } else if (contiene_T && contiene_U) {
    cout << "Invalida (contiene T y U)." << endl;
  } else if (solo_caracteres_ACGT && !contiene_U) { // Prioridad para clasificar como ADN
    cout << "ADN";
    // Si la longitud de la secuencia es múltiplo de 3, traducir todos los codones
    // if (!secuencia_procesada.empty() && secuencia_procesada.length() % 3 == 0) {
    //  string traduccion_codones_str;
    //  for (size_t i = 0; i < secuencia_procesada.length(); i += 3) {
    //    // Extraer el codón actual de 3 bases
    //    string codon = secuencia_procesada.substr(i, 3);
    //    traduccion_codones_str += obtenerAminoacidoDeCodon(codon, true); // true indica que es ADN
    //    if (i + 3 < secuencia_procesada.length()) {                      // Añadir coma si no es el último codón
    //      traduccion_codones_str += ", ";
    //    }
    //  }
    //  cout << " (Codones: " << traduccion_codones_str << ")";
    //}
    cout << endl;
  } else if (solo_caracteres_ACGU && !contiene_T) { // Siguiente prioridad: ARN
    cout << "ARN";
    // Si la longitud de la secuencia es múltiplo de 3, traducir todos los codones
    if (!secuencia_procesada.empty() && secuencia_procesada.length() % 3 == 0) {
      string traduccion_codones_str;
      for (size_t i = 0; i < secuencia_procesada.length(); i += 3) {
        // Extraer el codón actual de 3 bases
        string codon = secuencia_procesada.substr(i, 3);
        traduccion_codones_str += obtenerAminoacidoDeCodon(codon, false); // false indica que es ARN
        if (i + 3 < secuencia_procesada.length()) {                       // Añadir coma si no es el último codón
          traduccion_codones_str += ", ";
        }
      }
      cout << " (Codones: " << traduccion_codones_str << ")";
    }
    cout << endl;
  } else if (solo_caracteres_SLA_proteina && !contiene_U) { // Última prioridad: Proteína
    cout << "Proteina";
    if (!secuencia_procesada.empty()) {
      string lista_aminoacidos_str;
      for (size_t i = 0; i < secuencia_procesada.length(); ++i) {
        lista_aminoacidos_str += obtenerNombreCompletoAminoacido(secuencia_procesada[i]);
        if (i < secuencia_procesada.length() - 1) {
          lista_aminoacidos_str += ", ";
        }
      }
      if (secuencia_procesada.length() == 1) {
        cout << " (Aminoacido: " << lista_aminoacidos_str << ")";
      } else {
        cout << " (Aminoacidos: " << lista_aminoacidos_str << ")";
      }
    }
    cout << endl;
  } else {
    cout << "Tipo no determinado o mixta." << endl;
  }
}

// Función principal del programa
int main() {
  ifstream archivoEntrada("secuencias.txt");
  if (!archivoEntrada.is_open()) {
    cerr << "Error al abrir el archivo secuencias.txt" << endl;
    return 1;
  }

  string linea;
  while (getline(archivoEntrada, linea)) {
    if (!linea.empty() && linea.back() == '\r') {
      linea.pop_back();
    }
    procesarSecuencia(linea);
  }

  archivoEntrada.close();
  // procesarSecuencia("AUGGCCAUUGUAA"); // Ejemplo de secuencia de ARN
  return 0;
}
