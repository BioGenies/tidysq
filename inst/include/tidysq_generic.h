#ifndef RCPP_tidysq_H_GENERIC_
#define RCPP_tidysq_H_GENERIC_

namespace tidysq {
    template<typename RAW_CONTENT>
    RAW_CONTENT pack_raws_internal(RAW_CONTENT &unpacked,
                                   const unsigned short alph_size) {
        unsigned int in_len = unpacked.size();

        RAW_CONTENT ret((alph_size * in_len + 7) / 8);
        unsigned int out_byte = 0;

        int i = 0;

        if (alph_size == 2) {
            for (; i + 8 <= in_len; i += 8) {
                ret[out_byte] = (unpacked[i]) |
                                (unpacked[i + 1] << 2) |
                                (unpacked[i + 2] << 4) |
                                (unpacked[i + 3] << 6);
                ret[out_byte + 1] = (unpacked[i + 4]) |
                                    (unpacked[i + 5] << 2) |
                                    (unpacked[i + 6] << 4) |
                                    (unpacked[i + 7] << 6);
                out_byte += 2;
            }
            switch (in_len - i) {
                case 7:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2) |
                                    (unpacked[i + 2] << 4) |
                                    (unpacked[i + 3] << 6);
                    ret[out_byte + 1] = (unpacked[i + 4]) |
                                        (unpacked[i + 5] << 2) |
                                        (unpacked[i + 6] << 4);
                    break;
                case 6:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2) |
                                    (unpacked[i + 2] << 4) |
                                    (unpacked[i + 3] << 6);
                    ret[out_byte + 1] = (unpacked[i + 4]) |
                                        (unpacked[i + 5] << 2);
                    break;
                case 5:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2) |
                                    (unpacked[i + 2] << 4) |
                                    (unpacked[i + 3] << 6);
                    ret[out_byte + 1] = (unpacked[i + 4]);
                    break;
                case 4:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2) |
                                    (unpacked[i + 2] << 4) |
                                    (unpacked[i + 3] << 6);
                    break;
                case 3:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2) |
                                    (unpacked[i + 2] << 4);
                    break;
                case 2:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 2);
                    break;
                case 1:
                    ret[out_byte] = (unpacked[i]);
                    break;
            }
        } else if (alph_size == 3) {
            for (; i + 8 <= in_len; i += 8) {
                ret[out_byte] = (unpacked[i]) |
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
                    ret[out_byte] = (unpacked[i]) |
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
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 3) |
                                    (unpacked[i + 2] << 6);
                    ret[out_byte + 1] = (unpacked[i + 2] >> 2) |
                                        (unpacked[i + 3] << 1) |
                                        (unpacked[i + 4] << 4) |
                                        (unpacked[i + 5] << 7);
                    ret[out_byte + 2] = (unpacked[i + 5] >> 1);
                    break;
                case 5:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 3) |
                                    (unpacked[i + 2] << 6);
                    ret[out_byte + 1] = (unpacked[i + 2] >> 2) |
                                        (unpacked[i + 3] << 1) |
                                        (unpacked[i + 4] << 4);
                    break;
                case 4:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 3) |
                                    (unpacked[i + 2] << 6);
                    ret[out_byte + 1] = (unpacked[i + 2] >> 2) |
                                        (unpacked[i + 3] << 1);
                    break;
                case 3:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 3) |
                                    (unpacked[i + 2] << 6);
                    ret[out_byte + 1] = (unpacked[i + 2] >> 2);
                    break;
                case 2:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 3);
                    break;
                case 1:
                    ret[out_byte] = (unpacked[i]);
                    break;
            }
        } else if (alph_size == 4) {
            for (; i + 8 <= in_len; i += 8) {
                ret[out_byte] = (unpacked[i]) |
                                (unpacked[i + 1] << 4);
                ret[out_byte + 1] = (unpacked[i + 2]) |
                                    (unpacked[i + 3] << 4);
                ret[out_byte + 2] = (unpacked[i + 4]) |
                                    (unpacked[i + 5] << 4);
                ret[out_byte + 3] = (unpacked[i + 6]) |
                                    (unpacked[i + 7] << 4);
                out_byte += 4;
            }
            switch (in_len - i) {
                case 7:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    ret[out_byte + 1] = (unpacked[i + 2]) |
                                        (unpacked[i + 3] << 4);
                    ret[out_byte + 2] = (unpacked[i + 4]) |
                                        (unpacked[i + 5] << 4);
                    ret[out_byte + 3] = (unpacked[i + 6]);
                    break;
                case 6:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    ret[out_byte + 1] = (unpacked[i + 2]) |
                                        (unpacked[i + 3] << 4);
                    ret[out_byte + 2] = (unpacked[i + 4]) |
                                        (unpacked[i + 5] << 4);
                    break;
                case 5:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    ret[out_byte + 1] = (unpacked[i + 2]) |
                                        (unpacked[i + 3] << 4);
                    ret[out_byte + 2] = (unpacked[i + 4]);
                    break;
                case 4:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    ret[out_byte + 1] = (unpacked[i + 2]) |
                                        (unpacked[i + 3] << 4);
                    break;
                case 3:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    ret[out_byte + 1] = (unpacked[i + 2]);
                    break;
                case 2:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 4);
                    break;
                case 1:
                    ret[out_byte] = (unpacked[i]);
                    break;
            }
        } else if (alph_size == 5) {
            for (; i + 8 <= in_len; i += 8) {
                ret[out_byte] = (unpacked[i]) |
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
                    ret[out_byte] = (unpacked[i]) |
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
                    ret[out_byte] = (unpacked[i]) |
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
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 5);
                    ret[out_byte + 1] = (unpacked[i + 1] >> 3) |
                                        (unpacked[i + 2] << 2) |
                                        (unpacked[i + 3] << 7);
                    ret[out_byte + 2] = (unpacked[i + 3] >> 1) |
                                        (unpacked[i + 4] << 4);
                    ret[out_byte + 3] = (unpacked[i + 4] >> 4);
                    break;
                case 4:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 5);
                    ret[out_byte + 1] = (unpacked[i + 1] >> 3) |
                                        (unpacked[i + 2] << 2) |
                                        (unpacked[i + 3] << 7);
                    ret[out_byte + 2] = (unpacked[i + 3] >> 1);
                    break;
                case 3:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 5);
                    ret[out_byte + 1] = (unpacked[i + 1] >> 3) |
                                        (unpacked[i + 2] << 2);
                    break;
                case 2:
                    ret[out_byte] = (unpacked[i]) |
                                    (unpacked[i + 1] << 5);
                    ret[out_byte + 1] = (unpacked[i + 1] >> 3);
                    break;
                case 1:
                    ret[out_byte] = (unpacked[i]);
                    break;
            }
        }

        return ret;
    }

    template<typename RAW_CONTENT>
    unsigned int get_out_len_internal(RAW_CONTENT &packed,
                                      const unsigned int alph_size) {
        const unsigned int in_len = packed.size();
        unsigned short out_shift = 0;
        unsigned int out_len = 0;
        if (in_len == 0) {
            return out_len;
        }
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
                ((in_len % 5 == 4) && ((last & 62) == 0)) ||
                ((in_len % 5 == 0) && ((last & 248) == 0))) {
                out_shift = 1;
            }
            out_len = in_len * 8 / 5 - out_shift;
        }
        return out_len;
    }

    template<typename RAW_CONTENT_IN, typename RAW_CONTENT_OUT = RAW_CONTENT_IN>
    RAW_CONTENT_OUT unpack_raws_internal(RAW_CONTENT_IN &packed,
                                         const unsigned short alph_size) {
        const unsigned int in_len = packed.size();
        unsigned short out_len = get_out_len_internal<RAW_CONTENT_IN>(packed, alph_size);
        if (in_len == 0) {
            return RAW_CONTENT_OUT(0);
        }
        RAW_CONTENT_OUT ret(out_len);
        unsigned int in_byte = 0;

        int i = 0;
        if (alph_size == 2) {
            for (; i + 8 <= out_len; i += 8) {
                ret[i] = (packed[in_byte]) & 3;
                ret[i + 1] = (packed[in_byte] >> 2) & 3;
                ret[i + 2] = (packed[in_byte] >> 4) & 3;
                ret[i + 3] = (packed[in_byte] >> 6) & 3;
                ret[i + 4] = (packed[in_byte + 1]) & 3;
                ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                ret[i + 6] = (packed[in_byte + 1] >> 4) & 3;
                ret[i + 7] = (packed[in_byte + 1] >> 6) & 3;
                in_byte += 2;
            }
            switch (out_len - i) {
                case 7:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    ret[i + 2] = (packed[in_byte] >> 4) & 3;
                    ret[i + 3] = (packed[in_byte] >> 6) & 3;
                    ret[i + 4] = (packed[in_byte + 1]) & 3;
                    ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                    ret[i + 6] = (packed[in_byte + 1] >> 4) & 3;
                    break;
                case 6:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    ret[i + 2] = (packed[in_byte] >> 4) & 3;
                    ret[i + 3] = (packed[in_byte] >> 6) & 3;
                    ret[i + 4] = (packed[in_byte + 1]) & 3;
                    ret[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                    break;
                case 5:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    ret[i + 2] = (packed[in_byte] >> 4) & 3;
                    ret[i + 3] = (packed[in_byte] >> 6) & 3;
                    ret[i + 4] = (packed[in_byte + 1]) & 3;
                    break;
                case 4:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    ret[i + 2] = (packed[in_byte] >> 4) & 3;
                    ret[i + 3] = (packed[in_byte] >> 6) & 3;
                    break;
                case 3:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    ret[i + 2] = (packed[in_byte] >> 4) & 3;
                    break;
                case 2:
                    ret[i] = (packed[in_byte]) & 3;
                    ret[i + 1] = (packed[in_byte] >> 2) & 3;
                    break;
                case 1:
                    ret[i] = (packed[in_byte]) & 3;
                    break;
            }

        } else if (alph_size == 3) {
            for (; i + 8 <= out_len; i += 8) {
                ret[i] = (packed[in_byte]) & 7;
                ret[i + 1] = (packed[in_byte] >> 3) & 7;
                ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                             ((packed[in_byte + 1] << 2) & 7);
                ret[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                ret[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                ret[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                             ((packed[in_byte + 2] << 1) & 7);
                ret[i + 6] = (packed[in_byte + 2] >> 2) & 7;
                ret[i + 7] = (packed[in_byte + 2] >> 5) & 7;
                in_byte += 3;
            }
            switch (out_len - i) {
                case 7:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                 ((packed[in_byte + 1] << 2) & 7);
                    ret[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    ret[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    ret[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 7);
                    ret[i + 6] = (packed[in_byte + 2] >> 2) & 7;
                    break;
                case 6:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                 ((packed[in_byte + 1] << 2) & 7);
                    ret[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    ret[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    ret[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 7);
                    break;
                case 5:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                 ((packed[in_byte + 1] << 2) & 7);
                    ret[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    ret[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    break;
                case 4:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                 ((packed[in_byte + 1] << 2) & 7);
                    ret[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    break;
                case 3:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    ret[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                 ((packed[in_byte + 1] << 2) & 7);
                    break;
                case 2:
                    ret[i] = (packed[in_byte]) & 7;
                    ret[i + 1] = (packed[in_byte] >> 3) & 7;
                    break;
                case 1:
                    ret[i] = (packed[in_byte]) & 7;
                    break;
            }
        } else if (alph_size == 4) {
            for (; i + 8 <= out_len; i += 8) {
                ret[i] = (packed[in_byte]) & 15;
                ret[i + 1] = (packed[in_byte] >> 4) & 15;
                ret[i + 2] = (packed[in_byte + 1]) & 15;
                ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                ret[i + 4] = (packed[in_byte + 2]) & 15;
                ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                ret[i + 6] = (packed[in_byte + 3]) & 15;
                ret[i + 7] = (packed[in_byte + 3] >> 4) & 15;
                in_byte += 4;
            }
            switch (out_len - i) {
                case 7:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    ret[i + 2] = (packed[in_byte + 1]) & 15;
                    ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    ret[i + 4] = (packed[in_byte + 2]) & 15;
                    ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                    ret[i + 6] = (packed[in_byte + 3]) & 15;
                    break;
                case 6:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    ret[i + 2] = (packed[in_byte + 1]) & 15;
                    ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    ret[i + 4] = (packed[in_byte + 2]) & 15;
                    ret[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                    break;
                case 5:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    ret[i + 2] = (packed[in_byte + 1]) & 15;
                    ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    ret[i + 4] = (packed[in_byte + 2]) & 15;
                    break;
                case 4:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    ret[i + 2] = (packed[in_byte + 1]) & 15;
                    ret[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    break;
                case 3:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    ret[i + 2] = (packed[in_byte + 1]) & 15;
                    break;
                case 2:
                    ret[i] = (packed[in_byte]) & 15;
                    ret[i + 1] = (packed[in_byte] >> 4) & 15;
                    break;
                case 1:
                    ret[i] = (packed[in_byte]) & 15;
                    break;
            }
        } else if (alph_size == 5) {
            for (; i + 8 <= out_len; i += 8) {
                ret[i] = (packed[in_byte]) & 31;
                ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                             ((packed[in_byte + 1] << 3) & 31);
                ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                ret[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                             ((packed[in_byte + 2] << 1) & 31);
                ret[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                             ((packed[in_byte + 3] << 4) & 31);
                ret[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                ret[i + 6] = ((packed[in_byte + 3] >> 6) & 3) |
                             ((packed[in_byte + 4] << 2) & 31);
                ret[i + 7] = (packed[in_byte + 4] >> 3) & 31;
                in_byte += 5;
            }
            switch (out_len - i) {
                case 7:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    ret[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 31);
                    ret[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                 ((packed[in_byte + 3] << 4) & 31);
                    ret[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                    ret[i + 6] = ((packed[in_byte + 3] >> 6) & 3) |
                                 ((packed[in_byte + 4] << 2) & 31);
                    break;
                case 6:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    ret[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 31);
                    ret[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                 ((packed[in_byte + 3] << 4) & 31);
                    ret[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                    break;
                case 5:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    ret[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 31);
                    ret[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                 ((packed[in_byte + 3] << 4) & 31);
                    break;
                case 4:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    ret[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                 ((packed[in_byte + 2] << 1) & 31);
                    break;
                case 3:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    ret[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    break;
                case 2:
                    ret[i] = (packed[in_byte]) & 31;
                    ret[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                 ((packed[in_byte + 1] << 3) & 31);
                    break;
                case 1:
                    ret[i] = (packed[in_byte]) & 31;
                    break;
            }
        }

        return ret;
    }

    template<typename VECTOR>
    unsigned short get_alph_size_internal(VECTOR &alph) {
        return ceil(log2(alph.size() + 2));
    }
}

#endif // RCPP_tidysq_H_GEN_