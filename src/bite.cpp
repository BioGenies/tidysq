#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/internal/util.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(Rcpp::List x, Rcpp::IntegerVector indices) {
  const Sq<RCPP> sq = importFromR(x, "!");
  Sq<RCPP> ret(sq.length(), sq.alphabet());
  const AlphSize alph_size = sq.alphabet().alphabetSize();
  bool warning_called = false;
  Rcpp::StringVector NA_warning = R_NilValue;
  
  for (LenSq i = 0; i < sq.length(); ++i) {
    const Sequence<RCPP> sequence = sq[i];
    Sequence<RCPP> out_sequence(
        internal::calculatePackedLength(indices.length(), sq.alphabet()),
        indices.length()
    );
    for (LenSq j = 0; j < indices.length(); ++j) {
      LenSq index = indices[j] - 1;
      ElemPacked bitten_element = 0xff >> (8 - alph_size);

      if (index <= sequence.originalLength()) {
        LenSq seq_lowest_bit_index = alph_size * index;
        LenSq seq_highest_bit_index = seq_lowest_bit_index + alph_size - 1;
        LenSq seq_lowest_byte_index = seq_lowest_bit_index / 8;
        LenSq seq_highest_byte_index = seq_highest_bit_index / 8;
        unsigned short seq_lowest_bit_in_byte_index = seq_lowest_bit_index % 8;
        
        bitten_element = bitten_element &
          ((sequence[seq_lowest_byte_index] >> seq_lowest_bit_in_byte_index) |
          (sequence[seq_highest_byte_index] << (8 - seq_lowest_bit_in_byte_index)));
      } else if (!warning_called) {
        NA_warning = "some sequences are subsetted with index bigger than length - NA introduced";
        warning_called = true;
      }
      
      LenSq out_lowest_bit_index = alph_size * j;
      LenSq out_highest_bit_index = out_lowest_bit_index + alph_size - 1;
      LenSq out_lowest_byte_index = out_lowest_bit_index / 8;
      LenSq out_highest_byte_index = out_highest_bit_index / 8;
      unsigned short out_lowest_bit_in_byte_index = out_lowest_bit_index % 8;

      out_sequence[out_lowest_byte_index] = out_sequence[out_lowest_byte_index] |
        (bitten_element << out_lowest_bit_in_byte_index);
      if (out_highest_byte_index != out_lowest_byte_index) {
        out_sequence[out_highest_byte_index] = out_sequence[out_highest_byte_index] |
          (bitten_element >> (8 - out_lowest_bit_in_byte_index));
      }
    }
    ret[i] = out_sequence;
  }
  return Rcpp::List::create(Rcpp::Named("warning") = NA_warning,
                            Rcpp::Named("sq") = ret.exportToR());
}
