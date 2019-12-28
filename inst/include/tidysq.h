#ifndef RCPP_tidysq_H_GEN_
#define RCPP_tidysq_H_GEN_

#include "tidysq_RcppExports.h"

namespace tidysq {
  unsigned char* unpack_sq_to_char(Rcpp::RawVector packed, 
                                   const unsigned short alph_size) {
    Rcpp::RawVector unpacked = unpack_raws(packed, alph_size);
    unsigned char* chars = &unpacked[0];
    return chars;
  }
}

#endif // RCPP_tidysq_H_GEN_

