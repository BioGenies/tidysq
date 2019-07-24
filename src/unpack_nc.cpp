#include <Rcpp.h>

inline char retranslate_cnuc(char code, char na_char) {
  char ret;
  switch (code) {
  case 0: ret = 0; break;
  case 1: ret = 'A'; break;
  case 2: ret = 'C'; break;
  case 3: ret = 'G'; break;
  case 4: ret = 'T'; break;
  case 5: ret = 'U'; break;
  case 6: ret = '-'; break;
  case 7: ret = na_char;
  }
  return ret;
}


//' Unpack raw bytes of clean nucleotides sequence
//'
//' @param UNPACKED \code{raw} vector
//' @param ALPH_SIZE \code{integer}
// [[Rcpp::export]]
Rcpp::RawVector unpack_nc_cnuc(Rcpp::RawVector PACKED, const char na_char) {
  const unsigned int ALPH_SIZE = 3;
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  }
  unsigned char last = PACKED[in_len - 1];
  if (in_len % 3 == 1) {
    if ((last & 56) == 0) {
      out_shift = 1;
    } else out_shift = 0;
  } else if (in_len % 3 == 2) {
    if ((last & 126) == 0) {
      out_shift = 2;
    } else if ((last & 112) == 0) {
      out_shift = 1;
    } else out_shift = 0;
  } else if (in_len % 3 == 0) {
    if ((last & 252) == 0) {
      out_shift = 2; 
    } else if ((last & 224) == 0) {
      out_shift = 1;
    } else out_shift = 0;
  }
  
  const unsigned int out_len = in_len * 8 / 3 - out_shift;

  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  for (; i + 8 <= out_len; i += 8) {
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    ret[i + 3] = retranslate_cnuc((PACKED[in_byte + 1] >>  1) & 7, na_char);
    ret[i + 4] = retranslate_cnuc((PACKED[in_byte + 1] >>  4) & 7, na_char);
    ret[i + 5] = retranslate_cnuc((PACKED[in_byte + 1] >>  7) & 1, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 2] <<  1) & 7, na_char);
    ret[i + 6] = retranslate_cnuc((PACKED[in_byte + 2] >>  2) & 7, na_char);
    ret[i + 7] = retranslate_cnuc((PACKED[in_byte + 2] >>  5) & 7, na_char);
    in_byte += 3;
  }
  switch (out_len - i) {
  case 7: 
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    ret[i + 3] = retranslate_cnuc((PACKED[in_byte + 1] >>  1) & 7, na_char);
    ret[i + 4] = retranslate_cnuc((PACKED[in_byte + 1] >>  4) & 7, na_char);
    ret[i + 5] = retranslate_cnuc((PACKED[in_byte + 1] >>  7) & 1, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 2] <<  1) & 7, na_char);
    ret[i + 6] = retranslate_cnuc((PACKED[in_byte + 2] >>  2) & 7, na_char);
    break;
  case 6:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    ret[i + 3] = retranslate_cnuc((PACKED[in_byte + 1] >>  1) & 7, na_char);
    ret[i + 4] = retranslate_cnuc((PACKED[in_byte + 1] >>  4) & 7, na_char);
    ret[i + 5] = retranslate_cnuc((PACKED[in_byte + 1] >>  7) & 1, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 2] <<  1) & 7, na_char);
    break;
  case 5:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    ret[i + 3] = retranslate_cnuc((PACKED[in_byte + 1] >>  1) & 7, na_char);
    ret[i + 4] = retranslate_cnuc((PACKED[in_byte + 1] >>  4) & 7, na_char);
    break;
  case 4:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    ret[i + 3] = retranslate_cnuc((PACKED[in_byte + 1] >>  1) & 7, na_char);
    break;
  case 3:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    ret[i + 2] = retranslate_cnuc((PACKED[in_byte    ] >>  6) & 3, na_char) |
                 retranslate_cnuc((PACKED[in_byte + 1] <<  2) & 7, na_char);
    break;
  case 2:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    ret[i + 1] = retranslate_cnuc((PACKED[in_byte    ] >>  3) & 7, na_char);
    break;
  case 1:
    ret[i    ] = retranslate_cnuc((PACKED[in_byte    ]      ) & 7, na_char);
    break;
  }

  return ret;
}