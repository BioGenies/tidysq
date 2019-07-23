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

//' Pack raw bytes, but better
//' 
//' @param UNPACKED \code{raw} vector
//' @param ALPH_SIZE \code{integer}
// [[Rcpp::export]]
Rcpp::RawVector pack3(Rcpp::RawVector UNPACKED, 
                      const unsigned short ALPH_SIZE) {
  const unsigned int IN_LEN = UNPACKED.size();
  Rcpp::RawVector ret(ALPH_SIZE * IN_LEN  / 8);
  unsigned int out_byte = 0;
  if (ALPH_SIZE == 3) {
    for(int i = 0; i < IN_LEN; i += 8) {
      ret[out_byte    ] = (UNPACKED[i]         ) | (UNPACKED[i + 1] << 3) | (UNPACKED[i + 2] << 6);
      ret[out_byte + 1] = (UNPACKED[i + 2] >> 2) | (UNPACKED[i + 3] << 1) | (UNPACKED[i + 4] << 4) | (UNPACKED[i + 5] << 7);
      ret[out_byte + 2] = (UNPACKED[i + 5] >> 1) | (UNPACKED[i + 6] << 2) | (UNPACKED[i + 7] << 5);
      out_byte += 3;
    }
  }
  return ret;
}