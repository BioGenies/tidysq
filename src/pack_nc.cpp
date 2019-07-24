#include <Rcpp.h>
#define BYTE_SIZE 8
#define ZERO 0

inline char translate_cnuc(char letter) {
  char ret;
  switch (letter) {
  case 'a': case 'A': ret = 1; break;
  case 'c': case 'C': ret = 2; break;
  case 'g': case 'G': ret = 3; break;
  case 't': case 'T': ret = 4; break;
  case 'u': case 'U': ret = 5; break;
  case '-':           ret = 6; break;
  default:            ret = 7;
  }
  return ret;
}

//' Pack raw bytes, but better
//' 
//' @param UNPACKED \code{raw} vector
//' @param ALPH_SIZE \code{integer}
// [[Rcpp::export]]
Rcpp::RawVector pack_nc_cnuc(Rcpp::RawVector UNPACKED) {
  const unsigned int ALPH_SIZE = 3;
  unsigned int in_len = UNPACKED.size();
  
  Rcpp::RawVector ret((ALPH_SIZE * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  for (; i + 8 <= in_len; i += 8) {
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2) | 
                        (translate_cnuc(UNPACKED[i + 3]) << 1) | 
                        (translate_cnuc(UNPACKED[i + 4]) << 4) | 
                        (translate_cnuc(UNPACKED[i + 5]) << 7);
    ret[out_byte + 2] = (translate_cnuc(UNPACKED[i + 5]) >> 1) | 
                        (translate_cnuc(UNPACKED[i + 6]) << 2) | 
                        (translate_cnuc(UNPACKED[i + 7]) << 5);
    out_byte += 3;
  }
  switch (in_len - i) {
  case 7:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2) | 
                        (translate_cnuc(UNPACKED[i + 3]) << 1) | 
                        (translate_cnuc(UNPACKED[i + 4]) << 4) | 
                        (translate_cnuc(UNPACKED[i + 5]) << 7);
    ret[out_byte + 2] = (translate_cnuc(UNPACKED[i + 5]) >> 1) | 
                        (translate_cnuc(UNPACKED[i + 6]) << 2);
    break;
  case 6:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2) | 
                        (translate_cnuc(UNPACKED[i + 3]) << 1) | 
                        (translate_cnuc(UNPACKED[i + 4]) << 4) | 
                        (translate_cnuc(UNPACKED[i + 5]) << 7);
    ret[out_byte + 2] = (translate_cnuc(UNPACKED[i + 5]) >> 1);
    break;
  case 5:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2) | 
                        (translate_cnuc(UNPACKED[i + 3]) << 1) | 
                        (translate_cnuc(UNPACKED[i + 4]) << 4);
    break;
  case 4:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2) | 
                        (translate_cnuc(UNPACKED[i + 3]) << 1);
    break;
  case 3:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3) | 
                        (translate_cnuc(UNPACKED[i + 2]) << 6);
    ret[out_byte + 1] = (translate_cnuc(UNPACKED[i + 2]) >> 2); 
    break;
  case 2:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     ) | 
                        (translate_cnuc(UNPACKED[i + 1]) << 3);
    break;
  case 1:
    ret[out_byte    ] = (translate_cnuc(UNPACKED[i]    )     );
    break;
  }
  
  return ret;
}