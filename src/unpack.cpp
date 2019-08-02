#include <Rcpp.h>
#include <string>

unsigned short get_alph_size(Rcpp::CharacterVector alph);
Rcpp::CharacterVector match_raws(Rcpp::RawVector letters,
                                 Rcpp::CharacterVector alph,
                                 Rcpp::CharacterVector na_char);
// [[Rcpp::export]]
Rcpp::RawVector unpack_raws(Rcpp::RawVector packed, 
                            const unsigned short alph_size) {
  const unsigned int in_len = packed.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  } 
  unsigned int out_len;
  unsigned char last = packed[in_len - 1];
  if (alph_size == 2) {
    if ((last & 252) == 0) {
      out_shift = 3;
    } else if ((last & 240) == 0) {
      out_shift = 2;
    } else if ((last & 192) == 0) {
      out_shift = 1;
    }
    out_len = in_len * 4 - out_shift;
  } else if (alph_size == 3) {
    if (in_len % 3 == 1) {
      if ((last & 56) == 0) {
        out_shift = 1;
      } 
    } else if (in_len % 3 == 2) {
      if ((last & 126) == 0) {
        out_shift = 2;
      } else if ((last & 112) == 0) {
        out_shift = 1;
      } 
    } else if (in_len % 3 == 0) {
      if ((last & 252) == 0) {
        out_shift = 2; 
      } else if ((last & 224) == 0) {
        out_shift = 1;
      }
    }
    out_len = in_len * 8 / 3 - out_shift;
  } else if (alph_size == 4) {
    if ((last & 240) == 0) {
      out_shift = 1;
    }
    out_len = in_len * 2 - out_shift;
  } else if (alph_size == 5) {
    if (((in_len % 5 == 2) && ((last & 124) == 0)) ||
        ((in_len % 5 == 4) && ((last & 62 ) == 0)) ||
        ((in_len % 5 == 0) && ((last & 248) == 0))) {
      out_shift = 1;
    } 
    out_len = in_len * 8 / 5 - out_shift;
  }

  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  if (alph_size == 2) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      ret[i + 3] = (packed[in_byte    ] >> 6) & 3;
      ret[i + 4] = (packed[in_byte + 1]     ) & 3;
      ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
      ret[i + 6] = (packed[in_byte + 1] >> 4) & 3;
      ret[i + 7] = (packed[in_byte + 1] >> 6) & 3;
      in_byte += 2;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      ret[i + 3] = (packed[in_byte    ] >> 6) & 3;
      ret[i + 4] = (packed[in_byte + 1]     ) & 3;
      ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
      ret[i + 6] = (packed[in_byte + 1] >> 4) & 3;
      break;
    case 6:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      ret[i + 3] = (packed[in_byte    ] >> 6) & 3;
      ret[i + 4] = (packed[in_byte + 1]     ) & 3;
      ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
      break;
    case 5:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      ret[i + 3] = (packed[in_byte    ] >> 6) & 3;
      ret[i + 4] = (packed[in_byte + 1]     ) & 3;
      break;
    case 4:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      ret[i + 3] = (packed[in_byte    ] >> 6) & 3;
      break;
    case 3:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      ret[i + 2] = (packed[in_byte    ] >> 4) & 3;
      break;
    case 2:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      ret[i + 1] = (packed[in_byte    ] >> 2) & 3;
      break;
    case 1:
      ret[i    ] = (packed[in_byte    ]     ) & 3;
      break;
    }
    
  } else if (alph_size == 3) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      ret[i + 3] = (packed[in_byte + 1] >>  1) & 7 ;
      ret[i + 4] = (packed[in_byte + 1] >>  4) & 7 ;
      ret[i + 5] =((packed[in_byte + 1] >>  7) & 1) |
                  ((packed[in_byte + 2] <<  1) & 7);
      ret[i + 6] = (packed[in_byte + 2] >>  2) & 7 ;
      ret[i + 7] = (packed[in_byte + 2] >>  5) & 7 ;
      in_byte += 3;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      ret[i + 3] = (packed[in_byte + 1] >>  1) & 7 ;
      ret[i + 4] = (packed[in_byte + 1] >>  4) & 7 ;
      ret[i + 5] =((packed[in_byte + 1] >>  7) & 1) |
                  ((packed[in_byte + 2] <<  1) & 7);
      ret[i + 6] = (packed[in_byte + 2] >>  2) & 7 ;
      break;
    case 6:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      ret[i + 3] = (packed[in_byte + 1] >>  1) & 7 ;
      ret[i + 4] = (packed[in_byte + 1] >>  4) & 7 ;
      ret[i + 5] =((packed[in_byte + 1] >>  7) & 1) |
                  ((packed[in_byte + 2] <<  1) & 7);
      break;
    case 5:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      ret[i + 3] = (packed[in_byte + 1] >>  1) & 7 ;
      ret[i + 4] = (packed[in_byte + 1] >>  4) & 7 ;
      break;
    case 4:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      ret[i + 3] = (packed[in_byte + 1] >>  1) & 7 ;
      break;
    case 3:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      ret[i + 2] =((packed[in_byte    ] >>  6) & 3) |
                  ((packed[in_byte + 1] <<  2) & 7);
      break;
    case 2:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      ret[i + 1] = (packed[in_byte    ] >>  3) & 7 ;
      break;
    case 1:
      ret[i    ] = (packed[in_byte    ]      ) & 7 ;
      break;
    }
  } else if (alph_size == 4) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
      ret[i + 4] = (packed[in_byte + 2]     ) & 15;
      ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
      ret[i + 6] = (packed[in_byte + 3]     ) & 15;
      ret[i + 7] = (packed[in_byte + 3] >> 4) & 15;
      in_byte += 4;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
      ret[i + 4] = (packed[in_byte + 2]     ) & 15;
      ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
      ret[i + 6] = (packed[in_byte + 3]     ) & 15;
      break;
    case 6:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
      ret[i + 4] = (packed[in_byte + 2]     ) & 15;
      ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
      break;
    case 5:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
      ret[i + 4] = (packed[in_byte + 2]     ) & 15;
      break;
    case 4:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
      break;
    case 3:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      ret[i + 2] = (packed[in_byte + 1]     ) & 15;
      break;
    case 2:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      ret[i + 1] = (packed[in_byte    ] >> 4) & 15;
      break;
    case 1:
      ret[i    ] = (packed[in_byte    ]     ) & 15;
      break;
    }
  } else if (alph_size == 5) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      ret[i + 3] =((packed[in_byte + 1] >>  7) & 1 ) |
                  ((packed[in_byte + 2] <<  1) & 31);
      ret[i + 4] =((packed[in_byte + 2] >>  4) & 15) |
                  ((packed[in_byte + 3] <<  4) & 31);
      ret[i + 5] = (packed[in_byte + 3] >>  1) & 31 ;
      ret[i + 6] =((packed[in_byte + 3] >>  6) & 3 ) |
                  ((packed[in_byte + 4] <<  2) & 31);
      ret[i + 7] = (packed[in_byte + 4] >>  3) & 31 ;
      in_byte += 5;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      ret[i + 3] =((packed[in_byte + 1] >>  7) & 1 ) |
                  ((packed[in_byte + 2] <<  1) & 31);
      ret[i + 4] =((packed[in_byte + 2] >>  4) & 15) |
                  ((packed[in_byte + 3] <<  4) & 31);
      ret[i + 5] = (packed[in_byte + 3] >>  1) & 31 ;
      ret[i + 6] =((packed[in_byte + 3] >>  6) & 3 ) |
                  ((packed[in_byte + 4] <<  2) & 31);
      break;
    case 6:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      ret[i + 3] =((packed[in_byte + 1] >>  7) & 1 ) |
                  ((packed[in_byte + 2] <<  1) & 31);
      ret[i + 4] =((packed[in_byte + 2] >>  4) & 15) |
                  ((packed[in_byte + 3] <<  4) & 31);
      ret[i + 5] = (packed[in_byte + 3] >>  1) & 31 ;
      break;
    case 5:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      ret[i + 3] =((packed[in_byte + 1] >>  7) & 1 ) |
                  ((packed[in_byte + 2] <<  1) & 31);
      ret[i + 4] =((packed[in_byte + 2] >>  4) & 15) |
                  ((packed[in_byte + 3] <<  4) & 31);
      break;
    case 4:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      ret[i + 3] =((packed[in_byte + 1] >>  7) & 1 ) |
                  ((packed[in_byte + 2] <<  1) & 31);
      break;
    case 3:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      ret[i + 2] = (packed[in_byte + 1] >>  2) & 31 ;
      break;
    case 2:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      ret[i + 1] =((packed[in_byte    ] >>  5) & 7 ) |
                  ((packed[in_byte + 1] <<  3) & 31);
      break;
    case 1:
      ret[i    ] = (packed[in_byte    ]      ) & 31 ;
      break;
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::IntegerVector unpack_ints(Rcpp::RawVector packed,
                                const unsigned short alph_size) {
  if ((packed.size() == 0) or 
      ((packed.size() == 1) and 
         Rcpp::is_true(Rcpp::all(packed[0] == Rcpp::RawVector(1))))) return Rcpp::IntegerVector(0);
  Rcpp::RawVector unpacked = unpack_raws(packed, alph_size);
  Rcpp::IntegerVector ret(unpacked);
  return ret;
}

// [[Rcpp::export]]
Rcpp::CharacterVector unpack_chars(Rcpp::RawVector packed,
                                   Rcpp::CharacterVector alph,
                                   Rcpp::CharacterVector na_char) {
  if ((packed.size() == 0) or 
      ((packed.size() == 1) and 
         Rcpp::is_true(Rcpp::all(packed[0] == Rcpp::RawVector(1))))) return Rcpp::CharacterVector(0);
  Rcpp::RawVector unpacked = unpack_raws(packed, get_alph_size(alph));
  Rcpp::CharacterVector ret = match_raws(unpacked, alph, na_char);
  return ret;
}

// [[Rcpp::export]]
Rcpp::CharacterVector unpack_string(Rcpp::RawVector packed,
                              Rcpp::CharacterVector alph,
                              const char na_char) {
  Rcpp::CharacterVector unpacked = unpack_chars(packed, alph, na_char);
  std::string ret = Rcpp::collapse(unpacked);
  return ret;
}