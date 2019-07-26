#include <Rcpp.h>

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

inline char translate_nuc(char letter) {
  char ret;
  switch (letter) {
  case 'a': case 'A': ret = 1 ; break;
  case 'c': case 'C': ret = 2 ; break;
  case 'g': case 'G': ret = 3 ; break;
  case 't': case 'T': ret = 4 ; break;
  case 'u': case 'U': ret = 5 ; break;
  case 'w': case 'W': ret = 6 ; break;
  case 's': case 'S': ret = 7 ; break;
  case 'm': case 'M': ret = 8 ; break;
  case 'k': case 'K': ret = 9 ; break;
  case 'r': case 'R': ret = 10; break;
  case 'y': case 'Y': ret = 11; break;
  case 'b': case 'B': ret = 12; break;
  case 'd': case 'D': ret = 13; break;
  case 'h': case 'H': ret = 14; break;
  case 'v': case 'V': ret = 15; break;
  case 'n': case 'N': ret = 16; break;
  case '-':           ret = 17; break;
  default:            ret = 31;
  }
  return ret;
}

inline char translate_cami(char letter) {
  char ret;
  switch (letter) {
  case 'a': case 'A': ret = 1 ; break;
  case 'c': case 'C': ret = 2 ; break;
  case 'd': case 'D': ret = 3 ; break;
  case 'e': case 'E': ret = 4 ; break;
  case 'f': case 'F': ret = 5 ; break;
  case 'g': case 'G': ret = 6 ; break;
  case 'h': case 'H': ret = 7 ; break;
  case 'i': case 'I': ret = 8 ; break;
  case 'k': case 'K': ret = 9 ; break;
  case 'l': case 'L': ret = 10; break;
  case 'm': case 'M': ret = 11; break;
  case 'n': case 'N': ret = 12; break;
  case 'p': case 'P': ret = 13; break;
  case 'q': case 'Q': ret = 14; break;
  case 'r': case 'R': ret = 15; break;
  case 's': case 'S': ret = 16; break;
  case 't': case 'T': ret = 17; break;
  case 'v': case 'V': ret = 18; break;
  case 'w': case 'W': ret = 19; break;
  case 'y': case 'Y': ret = 20; break;
  case '-':           ret = 21; break;
  default:            ret = 31;
  }
  return ret;
}

inline char translate_ami(char letter) {
  char ret;
  switch (letter) {
  case 'a': case 'A': ret = 1 ; break;
  case 'b': case 'B': ret = 2 ; break;
  case 'c': case 'C': ret = 3 ; break;
  case 'd': case 'D': ret = 4 ; break;
  case 'e': case 'E': ret = 5 ; break;
  case 'f': case 'F': ret = 6 ; break;
  case 'g': case 'G': ret = 7 ; break;
  case 'h': case 'H': ret = 8 ; break;
  case 'i': case 'I': ret = 9 ; break;
  case 'j': case 'J': ret = 10; break;
  case 'k': case 'K': ret = 11; break;
  case 'l': case 'L': ret = 12; break;
  case 'm': case 'M': ret = 13; break;
  case 'n': case 'N': ret = 14; break;
  case 'o': case 'O': ret = 15; break;
  case 'p': case 'P': ret = 16; break;
  case 'q': case 'Q': ret = 17; break;
  case 'r': case 'R': ret = 18; break;
  case 's': case 'S': ret = 19; break;
  case 't': case 'T': ret = 20; break;
  case 'u': case 'U': ret = 21; break;
  case 'v': case 'V': ret = 22; break;
  case 'w': case 'W': ret = 23; break;
  case 'x': case 'X': ret = 24; break;
  case 'y': case 'Y': ret = 25; break;
  case 'z': case 'Z': ret = 26; break;
  case '-':           ret = 27; break;
  default:            ret = 31;
  }
  return ret;
}

inline char translate_unt(char letter, Rcpp::RawVector alph, const unsigned short ALPH_SIZE) {
  char ret = (1 << ALPH_SIZE) - 1;
  for (int i = 0; i < alph.size(); i ++) {
    if (alph[i] == letter) {
      ret = i + 1;
      break;
    }
  }
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector pack_cnuc(Rcpp::RawVector UNPACKED) {
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

// [[Rcpp::export]]
Rcpp::RawVector pack_nuc(Rcpp::RawVector UNPACKED) {
  const unsigned int ALPH_SIZE = 5;
  unsigned int in_len = UNPACKED.size();
  
  Rcpp::RawVector ret((ALPH_SIZE * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  for (; i + 8 <= in_len; i += 8) {
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2) | 
                        (translate_nuc(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_nuc(UNPACKED[i + 3]) >> 1) | 
                        (translate_nuc(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_nuc(UNPACKED[i + 4]) >> 4) |
                        (translate_nuc(UNPACKED[i + 5]) << 1) |
                        (translate_nuc(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_nuc(UNPACKED[i + 6]) >> 2) |
                        (translate_nuc(UNPACKED[i + 7]) << 3);
    out_byte += 5;
  }
  switch (in_len - i) {
  case 7:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2) | 
                        (translate_nuc(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_nuc(UNPACKED[i + 3]) >> 1) | 
                        (translate_nuc(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_nuc(UNPACKED[i + 4]) >> 4) |
                        (translate_nuc(UNPACKED[i + 5]) << 1) |
                        (translate_nuc(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_nuc(UNPACKED[i + 6]) >> 2);
    break;
  case 6:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2) | 
                        (translate_nuc(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_nuc(UNPACKED[i + 3]) >> 1) | 
                        (translate_nuc(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_nuc(UNPACKED[i + 4]) >> 4) |
                        (translate_nuc(UNPACKED[i + 5]) << 1);
    break;
  case 5:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2) | 
                        (translate_nuc(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_nuc(UNPACKED[i + 3]) >> 1) | 
                        (translate_nuc(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_nuc(UNPACKED[i + 4]) >> 4);
    break;
  case 4:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2) | 
                        (translate_nuc(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_nuc(UNPACKED[i + 3]) >> 1);
    break;
  case 3:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3) | 
                        (translate_nuc(UNPACKED[i + 2]) << 2);
    break;
  case 2:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     ) | 
                        (translate_nuc(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_nuc(UNPACKED[i + 1]) >> 3);
    break;
  case 1:
    ret[out_byte    ] = (translate_nuc(UNPACKED[i]    )     );
    break;
  }
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector pack_cami(Rcpp::RawVector UNPACKED) {
  const unsigned int ALPH_SIZE = 5;
  unsigned int in_len = UNPACKED.size();
  
  Rcpp::RawVector ret((ALPH_SIZE * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  for (; i + 8 <= in_len; i += 8) {
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2) | 
                        (translate_cami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_cami(UNPACKED[i + 3]) >> 1) | 
                        (translate_cami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_cami(UNPACKED[i + 4]) >> 4) |
                        (translate_cami(UNPACKED[i + 5]) << 1) |
                        (translate_cami(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_cami(UNPACKED[i + 6]) >> 2) |
                        (translate_cami(UNPACKED[i + 7]) << 3);
    out_byte += 5;
  }
  switch (in_len - i) {
  case 7:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2) | 
                        (translate_cami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_cami(UNPACKED[i + 3]) >> 1) | 
                        (translate_cami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_cami(UNPACKED[i + 4]) >> 4) |
                        (translate_cami(UNPACKED[i + 5]) << 1) |
                        (translate_cami(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_cami(UNPACKED[i + 6]) >> 2);
    break;
  case 6:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2) | 
                        (translate_cami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_cami(UNPACKED[i + 3]) >> 1) | 
                        (translate_cami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_cami(UNPACKED[i + 4]) >> 4) |
                        (translate_cami(UNPACKED[i + 5]) << 1);
    break;
  case 5:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2) | 
                        (translate_cami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_cami(UNPACKED[i + 3]) >> 1) | 
                        (translate_cami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_cami(UNPACKED[i + 4]) >> 4);
    break;
  case 4:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2) | 
                        (translate_cami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_cami(UNPACKED[i + 3]) >> 1);
    break;
  case 3:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3) | 
                        (translate_cami(UNPACKED[i + 2]) << 2);
    break;
  case 2:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     ) | 
                        (translate_cami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_cami(UNPACKED[i + 1]) >> 3);
      break;
  case 1:
    ret[out_byte    ] = (translate_cami(UNPACKED[i]    )     );
    break;
  }
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector pack_ami(Rcpp::RawVector UNPACKED) {
  const unsigned int ALPH_SIZE = 5;
  unsigned int in_len = UNPACKED.size();
  
  Rcpp::RawVector ret((ALPH_SIZE * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  for (; i + 8 <= in_len; i += 8) {
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2) | 
                        (translate_ami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_ami(UNPACKED[i + 3]) >> 1) | 
                        (translate_ami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_ami(UNPACKED[i + 4]) >> 4) |
                        (translate_ami(UNPACKED[i + 5]) << 1) |
                        (translate_ami(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_ami(UNPACKED[i + 6]) >> 2) |
                        (translate_ami(UNPACKED[i + 7]) << 3);
    out_byte += 5;
  }
  switch (in_len - i) {
  case 7:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2) | 
                        (translate_ami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_ami(UNPACKED[i + 3]) >> 1) | 
                        (translate_ami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_ami(UNPACKED[i + 4]) >> 4) |
                        (translate_ami(UNPACKED[i + 5]) << 1) |
                        (translate_ami(UNPACKED[i + 6]) << 6);
    ret[out_byte + 4] = (translate_ami(UNPACKED[i + 6]) >> 2);
    break;
  case 6:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2) | 
                        (translate_ami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_ami(UNPACKED[i + 3]) >> 1) | 
                        (translate_ami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_ami(UNPACKED[i + 4]) >> 4) |
                        (translate_ami(UNPACKED[i + 5]) << 1);
    break;
  case 5:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2) | 
                        (translate_ami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_ami(UNPACKED[i + 3]) >> 1) | 
                        (translate_ami(UNPACKED[i + 4]) << 4); 
    ret[out_byte + 3] = (translate_ami(UNPACKED[i + 4]) >> 4);
    break;
  case 4:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2) | 
                        (translate_ami(UNPACKED[i + 3]) << 7); 
    ret[out_byte + 2] = (translate_ami(UNPACKED[i + 3]) >> 1);
    break;
  case 3:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3) | 
                        (translate_ami(UNPACKED[i + 2]) << 2);
    break;
  case 2:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     ) | 
                        (translate_ami(UNPACKED[i + 1]) << 5); 
    ret[out_byte + 1] = (translate_ami(UNPACKED[i + 1]) >> 3);
      break;
  case 1:
    ret[out_byte    ] = (translate_ami(UNPACKED[i]    )     );
    break;
  }
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector pack_unt(Rcpp::RawVector UNPACKED, Rcpp::RawVector alph, 
                         const unsigned short ALPH_SIZE) {
  unsigned int in_len = UNPACKED.size();
  
  Rcpp::RawVector ret((ALPH_SIZE * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  
  if (ALPH_SIZE == 2) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 6);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 7], alph, ALPH_SIZE) << 6);
      out_byte += 2;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 6);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 4);
      break;
    case 6:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 6);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 2);
      break;
    case 5:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 6);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     );
      break;
    case 4:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 6);
      break;
    case 3:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 4); 
      break;
    case 2:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 2);
      break;
    case 1:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     );
      break;
    }
  } else if (ALPH_SIZE == 3) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 7);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 7], alph, ALPH_SIZE) << 5);
      out_byte += 3;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 7);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 2);
      break;
    case 6:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 7);
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) >> 1);
      break;
    case 5:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4);
      break;
    case 4:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 1);
      break;
    case 3:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 6);
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) >> 2); 
      break;
    case 2:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 3);
      break;
    case 1:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     );
      break;
    }
  } else if (ALPH_SIZE == 4) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 7], alph, ALPH_SIZE) << 4); 
      out_byte += 4;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE)     );
      break;
    case 6:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 4);
      break;
    case 5:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE)     );
      break;
    case 4:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 4);
      break;
    case 3:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE)     );
      break;
    case 2:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 4);
        break;
    case 1:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     );
      break;
    }
  } else if (ALPH_SIZE == 5) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 7); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) >> 4) |
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 1) |
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 6);
      ret[out_byte + 4] = (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) >> 2) |
                          (translate_unt(UNPACKED[i + 7], alph, ALPH_SIZE) << 3);
      out_byte += 5;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 7); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) >> 4) |
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 1) |
                          (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) << 6);
      ret[out_byte + 4] = (translate_unt(UNPACKED[i + 6], alph, ALPH_SIZE) >> 2);
      break;
    case 6:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 7); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) >> 4) |
                          (translate_unt(UNPACKED[i + 5], alph, ALPH_SIZE) << 1);
      break;
    case 5:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 7); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) >> 1) | 
                          (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) << 4); 
      ret[out_byte + 3] = (translate_unt(UNPACKED[i + 4], alph, ALPH_SIZE) >> 4);
      break;
    case 4:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2) | 
                          (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) << 7); 
      ret[out_byte + 2] = (translate_unt(UNPACKED[i + 3], alph, ALPH_SIZE) >> 1);
      break;
    case 3:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3) | 
                          (translate_unt(UNPACKED[i + 2], alph, ALPH_SIZE) << 2);
      break;
    case 2:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     ) | 
                          (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) << 5); 
      ret[out_byte + 1] = (translate_unt(UNPACKED[i + 1], alph, ALPH_SIZE) >> 3);
        break;
    case 1:
      ret[out_byte    ] = (translate_unt(UNPACKED[i]    , alph, ALPH_SIZE)     );
      break;
    }
  } 
  
  return ret;
}