#include <Rcpp.h>
#include <iostream>
#include <bitset>

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

inline char retranslate_nuc(char code, char na_char) {
  char ret;
  switch (code) {
  case 0 : ret = 0; break;
  case 1 : ret = 'A'; break;
  case 2 : ret = 'C'; break;
  case 3 : ret = 'G'; break;
  case 4 : ret = 'T'; break;
  case 5 : ret = 'U'; break;
  case 6 : ret = 'W'; break;
  case 7 : ret = 'S'; break;
  case 8 : ret = 'M'; break;
  case 9 : ret = 'K'; break;
  case 10: ret = 'R'; break; 
  case 11: ret = 'Y'; break;
  case 12: ret = 'B'; break;
  case 13: ret = 'D'; break;
  case 14: ret = 'H'; break;
  case 15: ret = 'V'; break;
  case 16: ret = 'N'; break;
  case 17: ret = '-'; break;
  default: ret = na_char;
  }
  return ret;
}

inline char retranslate_cami(char code, char na_char) {
  char ret;
  switch (code) {
  case 0 : ret = 0; break;
  case 1 : ret = 'A'; break;
  case 2 : ret = 'C'; break;
  case 3 : ret = 'D'; break;
  case 4 : ret = 'E'; break;
  case 5 : ret = 'F'; break;
  case 6 : ret = 'G'; break;
  case 7 : ret = 'H'; break;
  case 8 : ret = 'I'; break;
  case 9 : ret = 'K'; break;
  case 10: ret = 'L'; break; 
  case 11: ret = 'M'; break;
  case 12: ret = 'N'; break;
  case 13: ret = 'P'; break;
  case 14: ret = 'Q'; break;
  case 15: ret = 'R'; break;
  case 16: ret = 'S'; break;
  case 17: ret = 'T'; break;
  case 18: ret = 'V'; break;
  case 19: ret = 'W'; break;
  case 20: ret = 'Y'; break;
  case 21: ret = '-'; break;
  default: ret = na_char;
  }
  return ret;
}

inline char retranslate_ami(char code, char na_char) {
  char ret;
  switch (code) {
  case 0 : ret = 0; break;
  case 1 : ret = 'A'; break;
  case 2 : ret = 'B'; break;
  case 3 : ret = 'C'; break;
  case 4 : ret = 'D'; break;
  case 5 : ret = 'E'; break;
  case 6 : ret = 'F'; break;
  case 7 : ret = 'G'; break;
  case 8 : ret = 'H'; break;
  case 9 : ret = 'I'; break;
  case 10: ret = 'J'; break;
  case 11: ret = 'K'; break;
  case 12: ret = 'L'; break; 
  case 13: ret = 'M'; break;
  case 14: ret = 'N'; break;
  case 15: ret = 'O'; break;
  case 16: ret = 'P'; break;
  case 17: ret = 'Q'; break;
  case 18: ret = 'R'; break;
  case 19: ret = 'S'; break;
  case 20: ret = 'T'; break;
  case 21: ret = 'U'; break;
  case 22: ret = 'V'; break;
  case 23: ret = 'W'; break;
  case 24: ret = 'X'; break;
  case 25: ret = 'Y'; break;
  case 26: ret = 'Z'; break;
  case 27: ret = '-'; break;
  default: ret = na_char;
  }
  return ret;
}

inline char retranslate_unt(char code, char na_char, Rcpp::RawVector alph, 
                            const unsigned short ALPH_SIZE) {
  char ret;
  if (code == 0) {
    ret = 0;
  } else if (code == (1 << ALPH_SIZE) - 1) {
    ret = na_char;
  } else {
    ret = alph[code - 1];
  }
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector unpack_cnuc(Rcpp::RawVector PACKED, const char na_char) {
  const unsigned int ALPH_SIZE = 3;
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  }
  unsigned char last = PACKED[in_len - 1];
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
  
  const unsigned int out_len = in_len * 8 / 3 - out_shift;
  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  for (; i + 8 <= out_len; i += 8) {
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    ret[i + 3] = retranslate_cnuc( (PACKED[in_byte + 1] >>  1) & 7 , na_char);
    ret[i + 4] = retranslate_cnuc( (PACKED[in_byte + 1] >>  4) & 7 , na_char);
    ret[i + 5] = retranslate_cnuc(((PACKED[in_byte + 1] >>  7) & 1) |
                                  ((PACKED[in_byte + 2] <<  1) & 7), na_char);
    ret[i + 6] = retranslate_cnuc( (PACKED[in_byte + 2] >>  2) & 7 , na_char);
    ret[i + 7] = retranslate_cnuc( (PACKED[in_byte + 2] >>  5) & 7 , na_char);
    in_byte += 3;
  }
  switch (out_len - i) {
  case 7: 
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    ret[i + 3] = retranslate_cnuc( (PACKED[in_byte + 1] >>  1) & 7 , na_char);
    ret[i + 4] = retranslate_cnuc( (PACKED[in_byte + 1] >>  4) & 7 , na_char);
    ret[i + 5] = retranslate_cnuc(((PACKED[in_byte + 1] >>  7) & 1) |
                                  ((PACKED[in_byte + 2] <<  1) & 7), na_char);
    ret[i + 6] = retranslate_cnuc( (PACKED[in_byte + 2] >>  2) & 7 , na_char);
    break;
  case 6:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    ret[i + 3] = retranslate_cnuc( (PACKED[in_byte + 1] >>  1) & 7 , na_char);
    ret[i + 4] = retranslate_cnuc( (PACKED[in_byte + 1] >>  4) & 7 , na_char);
    ret[i + 5] = retranslate_cnuc(((PACKED[in_byte + 1] >>  7) & 1) |
                                  ((PACKED[in_byte + 2] <<  1) & 7), na_char);
    break;
  case 5:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    ret[i + 3] = retranslate_cnuc( (PACKED[in_byte + 1] >>  1) & 7 , na_char);
    ret[i + 4] = retranslate_cnuc( (PACKED[in_byte + 1] >>  4) & 7 , na_char);
    break;
  case 4:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    ret[i + 3] = retranslate_cnuc( (PACKED[in_byte + 1] >>  1) & 7 , na_char);
    break;
  case 3:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    ret[i + 2] = retranslate_cnuc(((PACKED[in_byte    ] >>  6) & 3) |
                                  ((PACKED[in_byte + 1] <<  2) & 7), na_char);
    break;
  case 2:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    ret[i + 1] = retranslate_cnuc( (PACKED[in_byte    ] >>  3) & 7 , na_char);
    break;
  case 1:
    ret[i    ] = retranslate_cnuc( (PACKED[in_byte    ]      ) & 7 , na_char);
    break;
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector unpack_nuc(Rcpp::RawVector PACKED, const char na_char) {
  const unsigned int ALPH_SIZE = 5;
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  }
  unsigned char last = PACKED[in_len - 1];
  if (((in_len % 5 == 2) && ((last & 124) == 0)) ||
      ((in_len % 5 == 4) && ((last & 62 ) == 0)) ||
      ((in_len % 5 == 0) && ((last & 249) == 0))) {
      out_shift = 1;
  } 
  
  const unsigned int out_len = in_len * 8 / 5 - out_shift;

  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  for (; i + 8 <= out_len; i += 8) {
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_nuc(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_nuc(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_nuc( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_nuc(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                 ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    ret[i + 7] = retranslate_nuc( (PACKED[in_byte + 4] >>  3) & 31 , na_char);
    in_byte += 5;
  }
  switch (out_len - i) {
  case 7: 
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_nuc(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_nuc(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_nuc( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_nuc(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                 ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    break;
  case 6:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_nuc(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_nuc(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_nuc( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    break;
  case 5:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_nuc(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_nuc(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    break;
  case 4:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_nuc(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    break;
  case 3:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_nuc( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    break;
  case 2:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_nuc(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    break;
  case 1:
    ret[i    ] = retranslate_nuc( (PACKED[in_byte    ]      ) & 31 , na_char);
    break;
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector unpack_cami(Rcpp::RawVector PACKED, const char na_char) {
  const unsigned int ALPH_SIZE = 5;
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  }
  unsigned char last = PACKED[in_len - 1];
  if (((in_len % 5 == 2) && ((last & 124) == 0)) ||
      ((in_len % 5 == 4) && ((last & 62 ) == 0)) ||
      ((in_len % 5 == 0) && ((last & 249) == 0))) {
      out_shift = 1;
  } 
  
  const unsigned int out_len = in_len * 8 / 5 - out_shift;

  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  for (; i + 8 <= out_len; i += 8) {
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_cami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                  ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_cami(((PACKED[in_byte + 2] >>  4) & 15) |
                                  ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_cami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_cami(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                  ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    ret[i + 7] = retranslate_cami( (PACKED[in_byte + 4] >>  3) & 31 , na_char);
    in_byte += 5;
  }
  switch (out_len - i) {
  case 7: 
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_cami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                  ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_cami(((PACKED[in_byte + 2] >>  4) & 15) |
                                  ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_cami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_cami(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                  ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    break;
  case 6:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_cami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                  ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_cami(((PACKED[in_byte + 2] >>  4) & 15) |
                                  ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_cami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    break;
  case 5:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_cami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                  ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_cami(((PACKED[in_byte + 2] >>  4) & 15) |
                                  ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    break;
  case 4:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_cami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                  ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    break;
  case 3:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_cami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    break;
  case 2:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_cami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                  ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    break;
  case 1:
    ret[i    ] = retranslate_cami( (PACKED[in_byte    ]      ) & 31 , na_char);
    break;
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector unpack_ami(Rcpp::RawVector PACKED, const char na_char) {
  const unsigned int ALPH_SIZE = 5;
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  }
  unsigned char last = PACKED[in_len - 1];
  if (((in_len % 5 == 2) && ((last & 124) == 0)) ||
      ((in_len % 5 == 4) && ((last & 62 ) == 0)) ||
      ((in_len % 5 == 0) && ((last & 249) == 0))) {
      out_shift = 1;
  } 
  
  const unsigned int out_len = in_len * 8 / 5 - out_shift;

  Rcpp::RawVector ret(out_len);
  unsigned int in_byte = 0;

  int i = 0;
  for (; i + 8 <= out_len; i += 8) {
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_ami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_ami(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_ami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_ami(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                 ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    ret[i + 7] = retranslate_ami( (PACKED[in_byte + 4] >>  3) & 31 , na_char);
    in_byte += 5;
  }
  switch (out_len - i) {
  case 7: 
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_ami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_ami(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_ami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    ret[i + 6] = retranslate_ami(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                 ((PACKED[in_byte + 4] <<  2) & 31), na_char);
    break;
  case 6:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_ami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_ami(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    ret[i + 5] = retranslate_ami( (PACKED[in_byte + 3] >>  1) & 31 , na_char);
    break;
  case 5:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_ami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    ret[i + 4] = retranslate_ami(((PACKED[in_byte + 2] >>  4) & 15) |
                                 ((PACKED[in_byte + 3] <<  4) & 31), na_char);
    break;
  case 4:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    ret[i + 3] = retranslate_ami(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                 ((PACKED[in_byte + 2] <<  1) & 31), na_char);
    break;
  case 3:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    ret[i + 2] = retranslate_ami( (PACKED[in_byte + 1] >>  2) & 31 , na_char);
    break;
  case 2:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    ret[i + 1] = retranslate_ami(((PACKED[in_byte    ] >>  5) & 7 ) |
                                 ((PACKED[in_byte + 1] <<  3) & 31), na_char);
    break;
  case 1:
    ret[i    ] = retranslate_ami( (PACKED[in_byte    ]      ) & 31 , na_char);
    break;
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector unpack_unt(Rcpp::RawVector PACKED, const char na_char,
                           Rcpp::RawVector alph, const unsigned short ALPH_SIZE) {
  const unsigned int in_len = PACKED.size();
  unsigned short out_shift = 0;
  if (in_len == 0) {
    return Rcpp::RawVector(0);
  } 
  unsigned int out_len;
  unsigned char last = PACKED[in_len - 1];
  if (ALPH_SIZE == 2) {
    if ((last & 252) == 0) {
      out_shift = 3;
    } else if ((last & 240) == 0) {
      out_shift = 2;
    } else if ((last & 92) == 0) {
      out_shift = 1;
    }
    out_len = in_len * 4 - out_shift;
  } else if (ALPH_SIZE == 3) {
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
  } else if (ALPH_SIZE == 4) {
    if ((last & 240) == 0) {
      out_shift = 1;
    }
    out_len = in_len * 2 - out_shift;
  } else if (ALPH_SIZE == 5) {
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
  if (ALPH_SIZE == 2) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte    ] >> 6) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 1]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 1] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 7] = retranslate_unt((PACKED[in_byte + 1] >> 6) & 3, na_char, alph, ALPH_SIZE);
      in_byte += 2;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte    ] >> 6) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 1]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 1] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 6:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte    ] >> 6) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 1]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 1] >> 2) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 5:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte    ] >> 6) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 1]     ) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 4:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte    ] >> 6) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 3:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte    ] >> 4) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 2:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 2) & 3, na_char, alph, ALPH_SIZE);
      break;
    case 1:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 3, na_char, alph, ALPH_SIZE);
      break;
    }
    
  } else if (ALPH_SIZE == 3) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt( (PACKED[in_byte + 1] >>  1) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt( (PACKED[in_byte + 1] >>  4) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1) |
                                   ((PACKED[in_byte + 2] <<  1) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt( (PACKED[in_byte + 2] >>  2) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 7] = retranslate_unt( (PACKED[in_byte + 2] >>  5) & 7 , na_char, alph, ALPH_SIZE);
      in_byte += 3;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt( (PACKED[in_byte + 1] >>  1) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt( (PACKED[in_byte + 1] >>  4) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1) |
                                   ((PACKED[in_byte + 2] <<  1) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt( (PACKED[in_byte + 2] >>  2) & 7 , na_char, alph, ALPH_SIZE);
      break;
    case 6:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt( (PACKED[in_byte + 1] >>  1) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt( (PACKED[in_byte + 1] >>  4) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1) |
                                   ((PACKED[in_byte + 2] <<  1) & 7), na_char, alph, ALPH_SIZE);
      break;
    case 5:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt( (PACKED[in_byte + 1] >>  1) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt( (PACKED[in_byte + 1] >>  4) & 7 , na_char, alph, ALPH_SIZE);
      break;
    case 4:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt( (PACKED[in_byte + 1] >>  1) & 7 , na_char, alph, ALPH_SIZE);
      break;
    case 3:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt(((PACKED[in_byte    ] >>  6) & 3) |
                                   ((PACKED[in_byte + 1] <<  2) & 7), na_char, alph, ALPH_SIZE);
      break;
    case 2:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt( (PACKED[in_byte    ] >>  3) & 7 , na_char, alph, ALPH_SIZE);
      break;
    case 1:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 7 , na_char, alph, ALPH_SIZE);
      break;
    }
  } else if (ALPH_SIZE == 4) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 2]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 2] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt((PACKED[in_byte + 3]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 7] = retranslate_unt((PACKED[in_byte + 3] >> 4) & 15, na_char, alph, ALPH_SIZE);
      in_byte += 4;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 2]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 2] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt((PACKED[in_byte + 3]     ) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 6:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 2]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt((PACKED[in_byte + 2] >> 4) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 5:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt((PACKED[in_byte + 2]     ) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 4:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt((PACKED[in_byte + 1] >> 4) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 3:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt((PACKED[in_byte + 1]     ) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 2:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt((PACKED[in_byte    ] >> 4) & 15, na_char, alph, ALPH_SIZE);
      break;
    case 1:
      ret[i    ] = retranslate_unt((PACKED[in_byte    ]     ) & 15, na_char, alph, ALPH_SIZE);
      break;
    }
  } else if (ALPH_SIZE == 5) {
    for (; i + 8 <= out_len; i += 8) {
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                   ((PACKED[in_byte + 2] <<  1) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt(((PACKED[in_byte + 2] >>  4) & 15) |
                                   ((PACKED[in_byte + 3] <<  4) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt( (PACKED[in_byte + 3] >>  1) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                   ((PACKED[in_byte + 4] <<  2) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 7] = retranslate_unt( (PACKED[in_byte + 4] >>  3) & 31 , na_char, alph, ALPH_SIZE);
      in_byte += 5;
    }
    switch (out_len - i) {
    case 7: 
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                   ((PACKED[in_byte + 2] <<  1) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt(((PACKED[in_byte + 2] >>  4) & 15) |
                                   ((PACKED[in_byte + 3] <<  4) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt( (PACKED[in_byte + 3] >>  1) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 6] = retranslate_unt(((PACKED[in_byte + 3] >>  6) & 3 ) |
                                   ((PACKED[in_byte + 4] <<  2) & 31), na_char, alph, ALPH_SIZE);
      break;
    case 6:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                   ((PACKED[in_byte + 2] <<  1) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt(((PACKED[in_byte + 2] >>  4) & 15) |
                                   ((PACKED[in_byte + 3] <<  4) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 5] = retranslate_unt( (PACKED[in_byte + 3] >>  1) & 31 , na_char, alph, ALPH_SIZE);
      break;
    case 5:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                   ((PACKED[in_byte + 2] <<  1) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 4] = retranslate_unt(((PACKED[in_byte + 2] >>  4) & 15) |
                                   ((PACKED[in_byte + 3] <<  4) & 31), na_char, alph, ALPH_SIZE);
      break;
    case 4:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 3] = retranslate_unt(((PACKED[in_byte + 1] >>  7) & 1 ) |
                                   ((PACKED[in_byte + 2] <<  1) & 31), na_char, alph, ALPH_SIZE);
      break;
    case 3:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      ret[i + 2] = retranslate_unt( (PACKED[in_byte + 1] >>  2) & 31 , na_char, alph, ALPH_SIZE);
      break;
    case 2:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      ret[i + 1] = retranslate_unt(((PACKED[in_byte    ] >>  5) & 7 ) |
                                   ((PACKED[in_byte + 1] <<  3) & 31), na_char, alph, ALPH_SIZE);
      break;
    case 1:
      ret[i    ] = retranslate_unt( (PACKED[in_byte    ]      ) & 31 , na_char, alph, ALPH_SIZE);
      break;
    }
  }
  
  return ret;
}