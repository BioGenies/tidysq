#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

unsigned short C_get_alph_size(Rcpp::CharacterVector alph);
Rcpp::RawVector C_unpack_raws(Rcpp::RawVector packed,
                              const unsigned short alph_size);
unsigned int C_get_out_len(Rcpp::RawVector packed,
                           const unsigned int alph_size);
void C_unpack_raws_safe(RcppParallel::RVector<unsigned char> packed,
                        RcppParallel::RVector<unsigned char> ret,
                        const unsigned short alph_size);

class UnpackingWorker: public RcppParallel::Worker {
public:
  std::vector<RcppParallel::RVector<unsigned char>> packed;
  std::vector<RcppParallel::RVector<unsigned char>> unpacked;
  unsigned int alphabetSize;
  
  UnpackingWorker(Rcpp::List sq) {
    packed = std::vector<RcppParallel::RVector<unsigned char>>();
    unpacked = std::vector<RcppParallel::RVector<unsigned char>>();
    alphabetSize = C_get_alph_size(sq.attr("alphabet"));

    for (int i = 0; i < sq.size(); i++) {
      Rcpp::RawVector v = sq[i];
      packed.push_back(RcppParallel::RVector<unsigned char>(v));
      unpacked.push_back(RcppParallel::RVector<unsigned char>(
          Rcpp::RawVector(C_get_out_len(Rcpp::RawVector(v), alphabetSize))));
    }
  }
  void operator()(std::size_t a, std::size_t b) {
    for(int i = a; i < b; ++i) {
      C_unpack_raws_safe(packed[i], unpacked[i], alphabetSize);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List C_unpack_sq_parallel(Rcpp::List sq) {
  UnpackingWorker worker(sq);
  RcppParallel::parallelFor(0, sq.size(), worker);
  for (int i = 0; i < sq.size(); i++) {
    sq[i] = Rcpp::RawVector(worker.unpacked[i].begin(), worker.unpacked[i].end());
  }
  sq.attr("class") = "list";
  return sq;
}