#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <string>
#include <tidysq_generic.h>

unsigned short C_get_alph_size(Rcpp::CharacterVector alph);
Rcpp::CharacterVector C_match(Rcpp::RawVector letters,
                              Rcpp::CharacterVector alph,
                              Rcpp::CharacterVector na_letter);

// [[Rcpp::export]]
Rcpp::RawVector C_unpack_raws(Rcpp::RawVector packed,
                              const unsigned short alph_size) {
  return tidysq::unpack_raws_internal<Rcpp::RawVector, Rcpp::RawVector>(packed, alph_size);  
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_unpack_ints(Rcpp::RawVector packed,
                                  const unsigned short alph_size) {
  if ((packed.size() == 0) or 
        ((packed.size() == 1) and 
           Rcpp::is_true(Rcpp::all(packed[0] == Rcpp::RawVector(1))))) return Rcpp::IntegerVector(0);
  Rcpp::RawVector unpacked = C_unpack_raws(packed, alph_size);
  Rcpp::IntegerVector ret(unpacked);
  return ret;
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_unpack_chars(Rcpp::RawVector packed,
                                     Rcpp::CharacterVector alph,
                                     Rcpp::CharacterVector na_letter) {
  if ((packed.size() == 0) or 
        ((packed.size() == 1) and 
           Rcpp::is_true(Rcpp::all(packed[0] == Rcpp::RawVector(1))))) return Rcpp::CharacterVector(0);
  Rcpp::RawVector unpacked = C_unpack_raws(packed, C_get_alph_size(alph));
  Rcpp::CharacterVector ret = C_match(unpacked, alph, na_letter);
  return ret;
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_unpack_string(Rcpp::RawVector packed,
                                      Rcpp::CharacterVector alph,
                                      Rcpp::CharacterVector na_letter) {
  Rcpp::CharacterVector unpacked = C_unpack_chars(packed, alph, na_letter);
  std::string ret = Rcpp::collapse(unpacked);
  return ret;
}

void C_unpack_raws_safe(RcppParallel::RVector<unsigned char> packed,
                        RcppParallel::RVector<unsigned char> ret,
                        const unsigned short alph_size) {
  const unsigned int in_len = packed.size();
  const unsigned int out_len = ret.size();
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
}