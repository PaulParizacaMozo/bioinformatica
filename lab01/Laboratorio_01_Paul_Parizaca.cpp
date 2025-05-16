#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
using namespace std;

// Constante global: Tabla de traducción de codones de ARN a nombres de aminoácidos.
// Un codón es una secuencia de tres nucleótidos.
const map<string, string> tablaCodonesARN = {
    {"UUU", "Fenilalanina"},    {"UUC", "Fenilalanina"},    {"UUA", "Leucina"},         {"UUG", "Leucina"},
    {"CUU", "Leucina"},         {"CUC", "Leucina"},         {"CUA", "Leucina"},         {"CUG", "Leucina"},
    {"AUU", "Isoleucina"},      {"AUC", "Isoleucina"},      {"AUA", "Isoleucina"},      {"AUG", "Metionina"},
    {"GUU", "Valina"},          {"GUC", "Valina"},          {"GUA", "Valina"},          {"GUG", "Valina"},
    {"UCU", "Serina"},          {"UCC", "Serina"},          {"UCA", "Serina"},          {"UCG", "Serina"},
    {"CCU", "Prolina"},         {"CCC", "Prolina"},         {"CCA", "Prolina"},         {"CCG", "Prolina"},
    {"ACU", "Treonina"},        {"ACC", "Treonina"},        {"ACA", "Treonina"},        {"ACG", "Treonina"},
    {"GCU", "Alanina"},         {"GCC", "Alanina"},         {"GCA", "Alanina"},         {"GCG", "Alanina"},
    {"UAU", "Tirosina"},        {"UAC", "Tirosina"},        {"UAA", "DETENCION"},       {"UAG", "DETENCION"},
    {"UGA", "DETENCION"},       {"CAU", "Histidina"},       {"CAC", "Histidina"},       {"CAA", "Glutamina"},
    {"CAG", "Glutamina"},       {"AAU", "Asparagina"},      {"AAC", "Asparagina"},      {"AAA", "Lisina"},
    {"AAG", "Lisina"},          {"GAU", "Acido Aspartico"}, {"GAC", "Acido Aspartico"}, {"GAA", "Acido Glutamico"},
    {"GAG", "Acido Glutamico"}, {"UGU", "Cisteina"},        {"UGC", "Cisteina"},        {"UGG", "Triptofano"},
    {"CGU", "Arginina"},        {"CGC", "Arginina"},        {"CGA", "Arginina"},        {"CGG", "Arginina"},
    {"AGU", "Serina"},          {"AGC", "Serina"},          {"AGA", "Arginina"},        {"AGG", "Arginina"},
    {"GGU", "Glicina"},         {"GGC", "Glicina"},         {"GGA", "Glicina"},         {"GGG", "Glicina"}};

// Mapa de Aminoacidos
const map<char, string> nombresAminoacidosSLA = {
    {'A', "Alanina"},   {'R', "Arginina"},        {'N', "Asparagina"}, {'D', "Acido Aspartico"}, {'C', "Cisteina"},
    {'Q', "Glutamina"}, {'E', "Acido Glutamico"}, {'G', "Glicina"},    {'H', "Histidina"},       {'I', "Isoleucina"},
    {'L', "Leucina"},   {'K', "Lisina"},          {'M', "Metionina"},  {'F', "Fenilalanina"},    {'P', "Prolina"},
    {'S', "Serina"},    {'T', "Treonina"},        {'V', "Valina"},     {'W', "Triptofano"},      {'Y', "Tirosina"}};

// Función: obtenerAminoacidoDeCodon
// Propósito: Traduce un codón (de ADN o ARN) al aminoácido correspondiente o señal de detención.

string obtenerAminoacidoDeCodon(const string &codon, bool esADN) {
  if (codon.length() != 3) {
    return "Longitud de codon invalida";
  }
  string codonARN = codon;

  auto iterador = tablaCodonesARN.find(codonARN);
  if (iterador != tablaCodonesARN.end()) {
    return iterador->second;
  }
  return "Codon desconocido";
}

// Función: obtenerNombreCompletoAminoacido
// Propósito: Obtiene el nombre completo de un aminoácido a partir de su código de una letra (SLA).

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

string procesarSecuencia(const string &secuencia_original) {
  string secuencia_procesada = secuencia_original;
  // Convertir toda la secuencia a mayúsculas para un análisis no sensible a mayúsculas/minúsculas
  transform(secuencia_procesada.begin(), secuencia_procesada.end(), secuencia_procesada.begin(), ::toupper);

  // string resultado = "Secuencia original: \"" + secuencia_original + "\" -> ";
  string resultado = "";

  if (secuencia_procesada.empty()) { // Comprobar si la secuencia está vacía
    return resultado + "Vacia.";
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
    return resultado + "Contiene caracteres no validos.";
  } else if (contiene_T && contiene_U) {
    return resultado + "Invalida (contiene T y U).";
  } else if (solo_caracteres_ACGT && !contiene_U) { // Prioridad para clasificar como ADN
    resultado += "ADN";
    return resultado;
  } else if (solo_caracteres_ACGU && !contiene_T) { // Siguiente prioridad: ARN
    resultado += "ARN";
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
      resultado += " (Codones: " + traduccion_codones_str + ")";
    }
    return resultado;
  } else if (solo_caracteres_SLA_proteina && !contiene_U) { // Última prioridad: Proteína
    resultado += "Proteina";
    if (!secuencia_procesada.empty()) {
      string lista_aminoacidos_str;
      for (size_t i = 0; i < secuencia_procesada.length(); ++i) {
        lista_aminoacidos_str += obtenerNombreCompletoAminoacido(secuencia_procesada[i]);
        if (i < secuencia_procesada.length() - 1) {
          lista_aminoacidos_str += ", ";
        }
      }
      if (secuencia_procesada.length() == 1) {
        resultado += " (Aminoacido: " + lista_aminoacidos_str + ")";
      } else {
        resultado += " (Aminoacidos: " + lista_aminoacidos_str + ")";
      }
    }
    return resultado;
  } else {
    return resultado + "Tipo no determinado o mixta.";
  }
}

// Pruebas para ADN
TEST_CASE("Pruebas de secuencias ADN") {

  // Caso básico con ADN sin errores
  SUBCASE("Secuencia ADN válida") { CHECK(procesarSecuencia("GATTACA") == "ADN"); }

  // Secuencia ADN válida pero con longitud múltiplo de 3
  SUBCASE("Secuencia ADN válida múltiplo de 3") { CHECK(procesarSecuencia("ACGT") == "ADN"); }

  // Secuencia ADN con un solo codón
  SUBCASE("Secuencia ADN con un codón válido") { CHECK(procesarSecuencia("GAT") == "ADN"); }

  // Secuencia ADN que no es divisible por 3 (longitud no divisible por 3)
  SUBCASE("Secuencia ADN no divisible por 3") {
    CHECK(procesarSecuencia("GATT") == "ADN");
    CHECK(procesarSecuencia("GATTA") == "ADN");
  }

  // Secuencia ADN válida con codones
  SUBCASE("Secuencia ADN con codones válidos") { CHECK(procesarSecuencia("ATGGTGCAATAG") == "ADN"); }

  // Secuencia ADN con caracteres no válidos
  SUBCASE("Secuencia ADN con caracteres no válidos") { CHECK(procesarSecuencia("ACGUX") == "Contiene caracteres no validos."); }

  // Secuencia ADN vacía
  SUBCASE("Secuencia ADN vacía") { CHECK(procesarSecuencia("") == "Vacia."); }

  // Secuencia ADN con codón de terminación
  SUBCASE("Secuencia ADN con codón de terminación") { CHECK(procesarSecuencia("TAA") == "ADN"); }

  // Secuencia ADN con una sola base
  SUBCASE("Secuencia ADN con una sola base") { CHECK(procesarSecuencia("A") == "ADN"); }

  // Secuencia ADN larga
  SUBCASE("Secuencia ADN larga") {
    CHECK(procesarSecuencia("GTAAACCCCTTTTCATTTAGACAGATCGACTCCTTATCCATTCTCAGAGATGTGTTGCTGGTCGCCG") == "ADN");
  }
}

// Pruebas para ARN
TEST_CASE("Pruebas de secuencias ARN") {

  // Secuencia ARN válida
  SUBCASE("Secuencia ARN válida") { CHECK(procesarSecuencia("AUGUCGUGA") == "ARN (Codones: Metionina, Serina, DETENCION)"); }

  // Secuencia ARN con codón de terminación
  SUBCASE("Secuencia ARN con codón de terminación") { CHECK(procesarSecuencia("UGA") == "ARN (Codones: DETENCION)"); }

  // Secuencia ARN con codones repetidos
  SUBCASE("Secuencia ARN con codones repetidos") { CHECK(procesarSecuencia("UUU") == "ARN (Codones: Fenilalanina)"); }

  // Secuencia ARN con longitud múltiplo de 3
  SUBCASE("Secuencia ARN válida múltiplo de 3") {
    CHECK(procesarSecuencia("AUGUGCCCGUAA") == "ARN (Codones: Metionina, Cisteina, Prolina, DETENCION)");
  }

  // Secuencia ARN con una letra (Proteína de una letra)
  SUBCASE("Secuencia ARN con una letra") { CHECK(procesarSecuencia("F") == "Proteina (Aminoacido: Fenilalanina)"); }

  // Secuencia ARN que no es divisible por 3
  SUBCASE("Secuencia ARN no divisible por 3") { CHECK(procesarSecuencia("UUUGA") == "ARN"); }

  // Secuencia ARN vacía
  SUBCASE("Secuencia ARN vacía") { CHECK(procesarSecuencia("") == "Vacia."); }

  // Secuencia ARN con caracteres no válidos
  SUBCASE("Secuencia ARN con caracteres no válidos") { CHECK(procesarSecuencia("GATO") == "Contiene caracteres no validos."); }

  // Secuencia ARN con codones de proteínas
  SUBCASE("Secuencia ARN codificando proteínas") {
    CHECK(procesarSecuencia("AUGUGCUGG") == "ARN (Codones: Metionina, Cisteina, Triptofano)");
  }

  // Secuencia ARN que contiene "T" (no válida en ARN)
  SUBCASE("Secuencia ARN con caracteres de ADN (contiene 'T')") {
    CHECK(procesarSecuencia("ACGTU") == "Invalida (contiene T y U).");
  }
}

// Pruebas para Proteinas
TEST_CASE("Pruebas de secuencias Proteínas") {

  // Secuencia Proteína de una sola letra (si no es ADN o ARN)
  SUBCASE("Secuencia Proteína de una sola letra") { CHECK(procesarSecuencia("F") == "Proteina (Aminoacido: Fenilalanina)"); }

  // Secuencia de proteínas con múltiples aminoácidos
  SUBCASE("Secuencia Proteína con varios aminoácidos") {
    CHECK(procesarSecuencia("RRR") == "Proteina (Aminoacidos: Arginina, Arginina, Arginina)");
  }

  // Secuencia de proteínas con caracteres válidos de aminoácidos
  SUBCASE("Secuencia Proteína con caracteres válidos") {
    CHECK(procesarSecuencia("HIV") == "Proteina (Aminoacidos: Histidina, Isoleucina, Valina)");
  }

  // Secuencia con caracteres no válidos (no pertenecen al mapa de aminoácidos)
  SUBCASE("Secuencia Proteína con caracteres no válidos") {
    CHECK(procesarSecuencia("X") == "Contiene caracteres no validos.");
  }

  // Secuencia Proteína vacía
  SUBCASE("Secuencia vacía para Proteína") { CHECK(procesarSecuencia("") == "Vacia."); }

  // Secuencia Proteína con una letra válida (si no es ARN o ADN)
  SUBCASE("Secuencia Proteína de una sola letra válida") {
    CHECK(procesarSecuencia("M") == "Proteina (Aminoacido: Metionina)");
  }

  // Secuencia ARN que debería ser procesada como ARN, no como proteína
  SUBCASE("Secuencia ARN que no debe ser procesada como Proteína") {
    CHECK(procesarSecuencia("AUG") == "ARN (Codones: Metionina)");
  }

  // Secuencia ADN que debería ser procesada como ADN, no como proteína
  SUBCASE("Secuencia ADN que no debe ser procesada como Proteína") { CHECK(procesarSecuencia("ACG") == "ADN"); }

  // Secuencia de proteína con longitud mayor a 3
  SUBCASE("Secuencia Proteína larga") {
    CHECK(procesarSecuencia("ARNDCEQ") ==
          "Proteina (Aminoacidos: Alanina, Arginina, Asparagina, Acido Aspartico, Cisteina, Acido Glutamico, Glutamina)");
  }
}

TEST_CASE("Secuencia ADN larga") {
  string secuenciaADN = "tgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttca"
                        "aggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg"
                        "caactccatcagagtggaaggaggcacctgggctgtgtatgaaaggcccaattttgctgg"
                        "gtacatgtacatcctaccccggggcgagtatcctgagtaccagcactggatgggcctcaa"
                        "cgaccgcctcagctcctgcagggctgttcacctgtctagtggaggccagtataagcttca"
                        "gatctttgagaaaggggattttaatggtcagatgcatgagaccacggaagactgcccttc"
                        "catcatggagcagttccacatgcgggaggtccactcctgtaaggtgctggagggcgcctg"
                        "gatcttctatgagctgcccaactaccgaggcaggcagtacctgctggacaagaaggagta"
                        "ccggaagcccgtcgactggggtgcagcttccccagctgtccagtctttccgccgcattgt"
                        "ggagtgatgatacagatgcggccaaacgctggctggccttgtcatccaaataagcattat"
                        "aaataaaacaattggcatgc";

  CHECK(procesarSecuencia(secuenciaADN) == "ADN");
}

TEST_CASE("Secuencias Proteínas largas") {

  // Cadena 1 de Proteína
  string secuenciaProteina1 = "MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFPTSTSLSPFYLR";
  string resultadoEsperado1 = "Proteina (Aminoacidos: Metionina, Acido Aspartico, Isoleucina, Alanina, "
                              "Isoleucina, Histidina, Histidina, Prolina, Triptofano, Isoleucina, Arginina, "
                              "Arginina, Prolina, Fenilalanina, Fenilalanina, Prolina, Fenilalanina, Histidina, "
                              "Serina, Prolina, Serina, Arginina, Leucina, Fenilalanina, Acido Aspartico, "
                              "Glutamina, Fenilalanina, Fenilalanina, Glicina, Acido Glutamico, Histidina, "
                              "Leucina, Leucina, Acido Glutamico, Serina, Acido Aspartico, Leucina, Fenilalanina, "
                              "Prolina, Treonina, Serina, Treonina, Serina, Leucina, Serina, Prolina, "
                              "Fenilalanina, Tirosina, Leucina, Arginina)";
  CHECK(procesarSecuencia(secuenciaProteina1) == resultadoEsperado1);

  // Cadena 2 de Proteína
  string secuenciaProteina2 = "PPSFLRAPSWIDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEV";
  string resultadoEsperado2 = "Proteina (Aminoacidos: Prolina, Prolina, Serina, Fenilalanina, Leucina, Arginina, "
                              "Alanina, Prolina, Serina, Triptofano, Isoleucina, Acido Aspartico, Treonina, "
                              "Glicina, Leucina, Serina, Acido Glutamico, Metionina, Arginina, Leucina, "
                              "Acido Glutamico, Lisina, Acido Aspartico, Arginina, Fenilalanina, Serina, "
                              "Valina, Asparagina, Leucina, Acido Aspartico, Valina, Lisina, Histidina, "
                              "Fenilalanina, Serina, Prolina, Acido Glutamico, Acido Glutamico, Leucina, "
                              "Lisina, Valina, Lisina, Valina, Leucina, Glicina, Acido Aspartico, Valina, "
                              "Isoleucina, Acido Glutamico, Valina)";
  CHECK(procesarSecuencia(secuenciaProteina2) == resultadoEsperado2);

  // Cadena 3 de Proteína
  string secuenciaProteina3 = "HGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRK";
  string resultadoEsperado3 = "Proteina (Aminoacidos: Histidina, Glicina, Lisina, Histidina, Acido Glutamico, "
                              "Acido Glutamico, Arginina, Glutamina, Acido Aspartico, Acido Glutamico, "
                              "Histidina, Glicina, Fenilalanina, Isoleucina, Serina, Arginina, Acido Glutamico, "
                              "Fenilalanina, Histidina, Arginina, Lisina, Tirosina, Arginina, Isoleucina, "
                              "Prolina, Alanina, Acido Aspartico, Valina, Acido Aspartico, Prolina, Leucina, "
                              "Treonina, Isoleucina, Treonina, Serina, Serina, Leucina, Serina, Serina, "
                              "Acido Aspartico, Glicina, Valina, Leucina, Treonina, Valina, Asparagina, Glicina, "
                              "Prolina, Arginina, Lisina)";
  CHECK(procesarSecuencia(secuenciaProteina3) == resultadoEsperado3);

  // Cadena 4 de Proteína
  string secuenciaProteina4 = "QAPGPERTIPITREEKPAVTAAPKK";
  string resultadoEsperado4 = "Proteina (Aminoacidos: Glutamina, Alanina, Prolina, Glicina, Prolina, Acido Glutamico, "
                              "Arginina, Treonina, Isoleucina, Prolina, Isoleucina, Treonina, Arginina, "
                              "Acido Glutamico, Acido Glutamico, Lisina, Prolina, Alanina, Valina, Treonina, "
                              "Alanina, Alanina, Prolina, Lisina, Lisina)";
  CHECK(procesarSecuencia(secuenciaProteina4) == resultadoEsperado4);
}
