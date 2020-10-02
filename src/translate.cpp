#include <Rcpp.h>
#include <map>

std::unordered_map<Rcpp::String, Rcpp::String> codonTable1 = {
  {"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "L"},
  {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},
  {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"}, {"ATG", "M"},
  {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},
  
  {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},
  {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},
  {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},
  {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},
  
  {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "*"}, {"TAG", "*"},
  {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
  {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"},
  {"GAT", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"},
  
  {"TGT", "C"}, {"TGC", "C"}, {"TGA", "*"}, {"TGG", "W"},
  {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},
  {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
  {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}
};

// [[Rcpp::export]]
Rcpp::CharacterVector extractCodons(std::string sequence) {
  Rcpp::CharacterVector ret(0);
  for (int i = 0; i <= (double)sequence.length() / 3 - 1; i++) {
    ret.push_back(sequence.substr(3 * i, 3));
  }
  return ret;
}

// [[Rcpp::export]]
Rcpp::String codonsToAminoAcids(Rcpp::CharacterVector codons) {
  Rcpp::String ret = "";
  for (Rcpp::CharacterVector::iterator i = codons.begin(); i != codons.end(); ++i) {
    ret += codonTable1[*i];
  }
  return ret;
}
