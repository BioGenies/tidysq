#include <Rcpp.h>
#define BYTE_SIZE 8
#define ZERO 0

//' Pack raw bytes
//' 
//' @param UNPACKED \code{raw} vector
//' @param ALPH_SIZE \code{integer}
// [[Rcpp::export]]
Rcpp::RawVector pack(Rcpp::RawVector UNPACKED, 
                     const unsigned short ALPH_SIZE) {
  const unsigned int IN_LEN = UNPACKED.size();
  Rcpp::RawVector ret((ALPH_SIZE * IN_LEN + BYTE_SIZE - 1) / BYTE_SIZE);
  unsigned int out_byte = ZERO;
  unsigned short bits_left = BYTE_SIZE;
  
  for (int i = ZERO; i < IN_LEN; i++) {
    if (bits_left >= ALPH_SIZE) {
      ret[out_byte] |= (UNPACKED[i] << (bits_left - ALPH_SIZE));
      bits_left -= ALPH_SIZE;
    } else {
      ret[out_byte] |= (UNPACKED[i] >> (ALPH_SIZE - bits_left));
      bits_left = ALPH_SIZE - bits_left;
      out_byte++;
      ret[out_byte] |= (UNPACKED[i] << (BYTE_SIZE - bits_left));
      bits_left = BYTE_SIZE - bits_left;
    }
  }
  return ret;
}

//' Unpack raw bytes
//' 
//' @param PACKED \code{raw} vector
//' @param ALPH_SIZE \code{integer}
// [[Rcpp::export]]
Rcpp::RawVector unpack(Rcpp::RawVector PACKED, 
                       const unsigned short ALPH_SIZE) {
  const unsigned int OUT_LEN = (PACKED.size() * BYTE_SIZE) / ALPH_SIZE;
  Rcpp::RawVector ret(OUT_LEN);
  const char MASK = (1u << ALPH_SIZE) - 1;
  unsigned int in_byte = ZERO;
  unsigned short bits_left = BYTE_SIZE;
  
  unsigned int i = ZERO;
  do {
    if (bits_left >= ALPH_SIZE) {
      ret[i] = ((PACKED[in_byte] >> (bits_left - ALPH_SIZE)) & MASK);
      if (ret[i] == 0) {
        for (int j = i; j < OUT_LEN; j++) {
          ret[j] = 0;
        }
        break;
      }
      bits_left -= ALPH_SIZE;
    } else {
      bits_left = ALPH_SIZE - bits_left;
      ret[i] = ((PACKED[in_byte] << bits_left) & MASK);
      in_byte++;
      ret[i] |= ((PACKED[in_byte] >> (BYTE_SIZE - bits_left)) & ((1u << bits_left) - 1));
      bits_left = BYTE_SIZE - bits_left;
    }
    i++;
  } while (i < OUT_LEN);
  return ret;
}
