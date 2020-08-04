#ifndef TIDYSQ_INTERNAL_UNPACK_H
#define TIDYSQ_INTERNAL_UNPACK_H

#include <stdexcept>

#include "types_general.h"
#include "internal_pack_util.h"

namespace tidysq {
    namespace internal {
        template<typename SEQUENCE_PACKED, 
                typename SEQUENCE_UNPACKED, 
                typename ALPHABET>
        inline SEQUENCE_UNPACKED unpack2(const SEQUENCE_PACKED& packed,
                                         const ALPHABET &alphabet) {
            lensq in_byte = 0;
            SEQUENCE_UNPACKED unpacked = reserveSpaceForUnpacked<SEQUENCE_PACKED, SEQUENCE_UNPACKED>(packed);
            lensq out_len = unpacked.size();

            lensq i = 0;
            for (; i + 8 <= out_len; i += 8) {
                unpacked[i] = (packed[in_byte]) & 3;
                unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                unpacked[i + 3] = (packed[in_byte] >> 6) & 3;
                unpacked[i + 4] = (packed[in_byte + 1]) & 3;
                unpacked[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                unpacked[i + 6] = (packed[in_byte + 1] >> 4) & 3;
                unpacked[i + 7] = (packed[in_byte + 1] >> 6) & 3;
                in_byte += 2;
            }
            switch (out_len - i) {
                case 7:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                    unpacked[i + 3] = (packed[in_byte] >> 6) & 3;
                    unpacked[i + 4] = (packed[in_byte + 1]) & 3;
                    unpacked[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                    unpacked[i + 6] = (packed[in_byte + 1] >> 4) & 3;
                    break;
                case 6:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                    unpacked[i + 3] = (packed[in_byte] >> 6) & 3;
                    unpacked[i + 4] = (packed[in_byte + 1]) & 3;
                    unpacked[i + 5] = (packed[in_byte + 1] >> 2) & 3;
                    break;
                case 5:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                    unpacked[i + 3] = (packed[in_byte] >> 6) & 3;
                    unpacked[i + 4] = (packed[in_byte + 1]) & 3;
                    break;
                case 4:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                    unpacked[i + 3] = (packed[in_byte] >> 6) & 3;
                    break;
                case 3:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    unpacked[i + 2] = (packed[in_byte] >> 4) & 3;
                    break;
                case 2:
                    unpacked[i] = (packed[in_byte]) & 3;
                    unpacked[i + 1] = (packed[in_byte] >> 2) & 3;
                    break;
                case 1:
                    unpacked[i] = (packed[in_byte]) & 3;
                    break;
            }
            return unpacked;
        }

        template<typename SEQUENCE_PACKED,
                typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        inline SEQUENCE_UNPACKED unpack3(const SEQUENCE_PACKED& packed,
                                         const ALPHABET &alphabet) {
            lensq in_byte = 0;
            SEQUENCE_UNPACKED unpacked = reserveSpaceForUnpacked<SEQUENCE_PACKED, SEQUENCE_UNPACKED>(packed);
            lensq out_len = unpacked.size();

            lensq i = 0;
            for (; i + 8 <= out_len; i += 8) {
                unpacked[i] = (packed[in_byte]) & 7;
                unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                  ((packed[in_byte + 1] << 2) & 7);
                unpacked[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                unpacked[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                unpacked[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                                  ((packed[in_byte + 2] << 1) & 7);
                unpacked[i + 6] = (packed[in_byte + 2] >> 2) & 7;
                unpacked[i + 7] = (packed[in_byte + 2] >> 5) & 7;
                in_byte += 3;
            }
            switch (out_len - i) {
                case 7:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                      ((packed[in_byte + 1] << 2) & 7);
                    unpacked[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    unpacked[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    unpacked[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 7);
                    unpacked[i + 6] = (packed[in_byte + 2] >> 2) & 7;
                    break;
                case 6:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                      ((packed[in_byte + 1] << 2) & 7);
                    unpacked[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    unpacked[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    unpacked[i + 5] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 7);
                    break;
                case 5:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                      ((packed[in_byte + 1] << 2) & 7);
                    unpacked[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    unpacked[i + 4] = (packed[in_byte + 1] >> 4) & 7;
                    break;
                case 4:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                      ((packed[in_byte + 1] << 2) & 7);
                    unpacked[i + 3] = (packed[in_byte + 1] >> 1) & 7;
                    break;
                case 3:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    unpacked[i + 2] = ((packed[in_byte] >> 6) & 3) |
                                      ((packed[in_byte + 1] << 2) & 7);
                    break;
                case 2:
                    unpacked[i] = (packed[in_byte]) & 7;
                    unpacked[i + 1] = (packed[in_byte] >> 3) & 7;
                    break;
                case 1:
                    unpacked[i] = (packed[in_byte]) & 7;
                    break;
            }
            return unpacked;
        }

        template<typename SEQUENCE_PACKED,
                typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        inline SEQUENCE_UNPACKED unpack4(const SEQUENCE_PACKED& packed,
                                         const ALPHABET &alphabet) {
            lensq in_byte = 0;
            SEQUENCE_UNPACKED unpacked = reserveSpaceForUnpacked<SEQUENCE_PACKED, SEQUENCE_UNPACKED>(packed);
            lensq out_len = unpacked.size();
            
            lensq i = 0;
            for (; i + 8 <= out_len; i += 8) {
                unpacked[i] = (packed[in_byte]) & 15;
                unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                unpacked[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                unpacked[i + 4] = (packed[in_byte + 2]) & 15;
                unpacked[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                unpacked[i + 6] = (packed[in_byte + 3]) & 15;
                unpacked[i + 7] = (packed[in_byte + 3] >> 4) & 15;
                in_byte += 4;
            }
            switch (out_len - i) {
                case 7:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                    unpacked[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    unpacked[i + 4] = (packed[in_byte + 2]) & 15;
                    unpacked[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                    unpacked[i + 6] = (packed[in_byte + 3]) & 15;
                    break;
                case 6:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                    unpacked[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    unpacked[i + 4] = (packed[in_byte + 2]) & 15;
                    unpacked[i + 5] = (packed[in_byte + 2] >> 4) & 15;
                    break;
                case 5:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                    unpacked[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    unpacked[i + 4] = (packed[in_byte + 2]) & 15;
                    break;
                case 4:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                    unpacked[i + 3] = (packed[in_byte + 1] >> 4) & 15;
                    break;
                case 3:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    unpacked[i + 2] = (packed[in_byte + 1]) & 15;
                    break;
                case 2:
                    unpacked[i] = (packed[in_byte]) & 15;
                    unpacked[i + 1] = (packed[in_byte] >> 4) & 15;
                    break;
                case 1:
                    unpacked[i] = (packed[in_byte]) & 15;
                    break;
            }
            return unpacked;
        }

        template<typename SEQUENCE_PACKED,
                typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        inline SEQUENCE_UNPACKED unpack5(const SEQUENCE_PACKED& packed,
                                         const ALPHABET &alphabet) {
            lensq in_byte = 0;
            SEQUENCE_UNPACKED unpacked = reserveSpaceForUnpacked<SEQUENCE_PACKED, SEQUENCE_UNPACKED>(packed);
            lensq out_len = unpacked.size();

            lensq i = 0;
            for (; i + 8 <= out_len; i += 8) {
                unpacked[i] = (packed[in_byte]) & 31;
                unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                  ((packed[in_byte + 1] << 3) & 31);
                unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                unpacked[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                  ((packed[in_byte + 2] << 1) & 31);
                unpacked[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                  ((packed[in_byte + 3] << 4) & 31);
                unpacked[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                unpacked[i + 6] = ((packed[in_byte + 3] >> 6) & 3) |
                                  ((packed[in_byte + 4] << 2) & 31);
                unpacked[i + 7] = (packed[in_byte + 4] >> 3) & 31;
                in_byte += 5;
            }
            switch (out_len - i) {
                case 7:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    unpacked[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 31);
                    unpacked[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                      ((packed[in_byte + 3] << 4) & 31);
                    unpacked[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                    unpacked[i + 6] = ((packed[in_byte + 3] >> 6) & 3) |
                                      ((packed[in_byte + 4] << 2) & 31);
                    break;
                case 6:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    unpacked[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 31);
                    unpacked[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                      ((packed[in_byte + 3] << 4) & 31);
                    unpacked[i + 5] = (packed[in_byte + 3] >> 1) & 31;
                    break;
                case 5:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    unpacked[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 31);
                    unpacked[i + 4] = ((packed[in_byte + 2] >> 4) & 15) |
                                      ((packed[in_byte + 3] << 4) & 31);
                    break;
                case 4:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    unpacked[i + 3] = ((packed[in_byte + 1] >> 7) & 1) |
                                      ((packed[in_byte + 2] << 1) & 31);
                    break;
                case 3:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    unpacked[i + 2] = (packed[in_byte + 1] >> 2) & 31;
                    break;
                case 2:
                    unpacked[i] = (packed[in_byte]) & 31;
                    unpacked[i + 1] = ((packed[in_byte] >> 5) & 7) |
                                      ((packed[in_byte + 1] << 3) & 31);
                    break;
                case 1:
                    unpacked[i] = (packed[in_byte]) & 31;
                    break;
            }
            return unpacked;
        }

        template<typename SEQUENCE_PACKED,
                typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        inline SEQUENCE_UNPACKED unpack(const SEQUENCE_PACKED& packed,
                                        const ALPHABET &alphabet) {
            switch (getAlphabetSize(alphabet)) {
                case 2:
                    return unpack2<SEQUENCE_PACKED, SEQUENCE_UNPACKED, ALPHABET>(packed, alphabet);
                case 3:
                    return unpack3<SEQUENCE_PACKED, SEQUENCE_UNPACKED, ALPHABET>(packed, alphabet);
                case 4:
                    return unpack4<SEQUENCE_PACKED, SEQUENCE_UNPACKED, ALPHABET>(packed, alphabet);
                case 5:
                    return unpack5<SEQUENCE_PACKED, SEQUENCE_UNPACKED, ALPHABET>(packed, alphabet);
                default:
                    throw std::invalid_argument("packed");
            }
        }
    }
}

#endif //TIDYSQ_INTERNAL_UNPACK_H
