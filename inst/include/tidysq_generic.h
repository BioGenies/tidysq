#ifndef RCPP_tidysq_H_GENERIC_
#define RCPP_tidysq_H_GENERIC_

namespace tidysq {
template<typename RAW_CONTENT>
RAW_CONTENT pack_raws_internal(RAW_CONTENT unpacked,
                               const unsigned short alph_size) {
  unsigned int in_len = unpacked.size();
  
  RAW_CONTENT ret((alph_size * in_len  + 7) / 8);
  unsigned int out_byte = 0;
  
  int i = 0;
  
  if (alph_size == 2) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4) | 
        (unpacked[i + 3] << 6);
      ret[out_byte + 1] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 2) | 
        (unpacked[i + 6] << 4) | 
        (unpacked[i + 7] << 6);
      out_byte += 2;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4) | 
        (unpacked[i + 3] << 6);
      ret[out_byte + 1] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 2) | 
        (unpacked[i + 6] << 4);
      break;
    case 6:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4) | 
        (unpacked[i + 3] << 6);
      ret[out_byte + 1] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 2);
      break;
    case 5:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4) | 
        (unpacked[i + 3] << 6);
      ret[out_byte + 1] = (unpacked[i + 4]     );
      break;
    case 4:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4) | 
        (unpacked[i + 3] << 6);
      break;
    case 3:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2) | 
        (unpacked[i + 2] << 4); 
      break;
    case 2:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 2);
      break;
    case 1:
      ret[out_byte    ] = (unpacked[i]         );
      break;
    }
  } else if (alph_size == 3) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2) | 
        (unpacked[i + 3] << 1) | 
        (unpacked[i + 4] << 4) | 
        (unpacked[i + 5] << 7);
      ret[out_byte + 2] = (unpacked[i + 5] >> 1) | 
        (unpacked[i + 6] << 2) | 
        (unpacked[i + 7] << 5);
      out_byte += 3;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2) | 
        (unpacked[i + 3] << 1) | 
        (unpacked[i + 4] << 4) | 
        (unpacked[i + 5] << 7);
      ret[out_byte + 2] = (unpacked[i + 5] >> 1) | 
        (unpacked[i + 6] << 2);
      break;
    case 6:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2) | 
        (unpacked[i + 3] << 1) | 
        (unpacked[i + 4] << 4) | 
        (unpacked[i + 5] << 7);
      ret[out_byte + 2] = (unpacked[i + 5] >> 1);
      break;
    case 5:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2) | 
        (unpacked[i + 3] << 1) | 
        (unpacked[i + 4] << 4);
      break;
    case 4:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2) | 
        (unpacked[i + 3] << 1);
      break;
    case 3:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3) | 
        (unpacked[i + 2] << 6);
      ret[out_byte + 1] = (unpacked[i + 2] >> 2); 
      break;
    case 2:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 3);
      break;
    case 1:
      ret[out_byte    ] = (unpacked[i]         );
      break;
    }
  } else if (alph_size == 4) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     ) | 
        (unpacked[i + 3] << 4); 
      ret[out_byte + 2] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 4); 
      ret[out_byte + 3] = (unpacked[i + 6]     ) | 
        (unpacked[i + 7] << 4); 
      out_byte += 4;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     ) | 
        (unpacked[i + 3] << 4); 
      ret[out_byte + 2] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 4); 
      ret[out_byte + 3] = (unpacked[i + 6]     );
      break;
    case 6:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     ) | 
        (unpacked[i + 3] << 4); 
      ret[out_byte + 2] = (unpacked[i + 4]     ) | 
        (unpacked[i + 5] << 4);
      break;
    case 5:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     ) | 
        (unpacked[i + 3] << 4); 
      ret[out_byte + 2] = (unpacked[i + 4]     );
      break;
    case 4:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     ) | 
        (unpacked[i + 3] << 4);
      break;
    case 3:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4); 
      ret[out_byte + 1] = (unpacked[i + 2]     );
      break;
    case 2:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 4);
      break;
    case 1:
      ret[out_byte    ] = (unpacked[i]         );
      break;
    }
  } else if (alph_size == 5) {
    for (; i + 8 <= in_len; i += 8) {
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2) | 
        (unpacked[i + 3] << 7); 
      ret[out_byte + 2] = (unpacked[i + 3] >> 1) | 
        (unpacked[i + 4] << 4); 
      ret[out_byte + 3] = (unpacked[i + 4] >> 4) |
        (unpacked[i + 5] << 1) |
        (unpacked[i + 6] << 6);
      ret[out_byte + 4] = (unpacked[i + 6] >> 2) |
        (unpacked[i + 7] << 3);
      out_byte += 5;
    }
    switch (in_len - i) {
    case 7:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2) | 
        (unpacked[i + 3] << 7); 
      ret[out_byte + 2] = (unpacked[i + 3] >> 1) | 
        (unpacked[i + 4] << 4); 
      ret[out_byte + 3] = (unpacked[i + 4] >> 4) |
        (unpacked[i + 5] << 1) |
        (unpacked[i + 6] << 6);
      ret[out_byte + 4] = (unpacked[i + 6] >> 2);
      break;
    case 6:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2) | 
        (unpacked[i + 3] << 7); 
      ret[out_byte + 2] = (unpacked[i + 3] >> 1) | 
        (unpacked[i + 4] << 4); 
      ret[out_byte + 3] = (unpacked[i + 4] >> 4) |
        (unpacked[i + 5] << 1);
      break;
    case 5:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2) | 
        (unpacked[i + 3] << 7); 
      ret[out_byte + 2] = (unpacked[i + 3] >> 1) | 
        (unpacked[i + 4] << 4); 
      ret[out_byte + 3] = (unpacked[i + 4] >> 4);
      break;
    case 4:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2) | 
        (unpacked[i + 3] << 7); 
      ret[out_byte + 2] = (unpacked[i + 3] >> 1);
      break;
    case 3:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3) | 
        (unpacked[i + 2] << 2);
      break;
    case 2:
      ret[out_byte    ] = (unpacked[i]         ) | 
        (unpacked[i + 1] << 5); 
      ret[out_byte + 1] = (unpacked[i + 1] >> 3);
      break;
    case 1:
      ret[out_byte    ] = (unpacked[i]         );
      break;
    }
  } 
  
  return ret;
}

template<typename VECTOR>
unsigned short get_alph_size_internal(VECTOR alph) {
  return ceil(log2(alph.size() + 2));
}
}

#endif // RCPP_tidysq_H_GEN_