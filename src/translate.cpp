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

std::unordered_map<int, std::unordered_map<Rcpp::String, Rcpp::String> > codonTables = {
  {1, codonTable1}
};

// [[Rcpp::export]]
Rcpp::CharacterVector Cpp_translate(std::vector<std::string> sq, int table) {
  Rcpp::CharacterVector ret(sq.size());
  for (std::size_t sq_index = 0; sq_index < sq.size(); sq_index++) {
    std::string sequence = sq[sq_index];
    for (int i = 0; i <= (double)sequence.length() / 3 - 1; i++) {
      ret[sq_index] += codonTables[table][sequence.substr(3 * i, 3)];
    }
  }
  return ret;
}
