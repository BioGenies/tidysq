#ifndef TIDYSQ_INTERNAL_PACK_H
#define TIDYSQ_INTERNAL_PACK_H

#include <stdexcept>

#include "types/general.h"
#include "internal_pack_util.h"

namespace tidysq {
    namespace internal {
        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED pack2(const SEQUENCE_UNPACKED& unpacked,
                                     const ALPHABET& alphabet) {
            lensq outByte = 0;
            lensq i = 0;
            SEQUENCE_PACKED packed = reserveSpaceForPacked<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
            for (; i + 8 <= unpacked.size(); i += 8) {
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4) |
                                  (unpacked[i + 3] << 6);
                packed[outByte + 1] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 2) |
                                      (unpacked[i + 6] << 4) |
                                      (unpacked[i + 7] << 6);
                outByte += 2;
            }
            switch (unpacked.size() - i) {
                case 7:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2) |
                                      (unpacked[i + 2] << 4) |
                                      (unpacked[i + 3] << 6);
                    packed[outByte + 1] = (unpacked[i + 4]) |
                                          (unpacked[i + 5] << 2) |
                                          (unpacked[i + 6] << 4);
                    break;
                case 6:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2) |
                                      (unpacked[i + 2] << 4) |
                                      (unpacked[i + 3] << 6);
                    packed[outByte + 1] = (unpacked[i + 4]) |
                                          (unpacked[i + 5] << 2);
                    break;
                case 5:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2) |
                                      (unpacked[i + 2] << 4) |
                                      (unpacked[i + 3] << 6);
                    packed[outByte + 1] = (unpacked[i + 4]);
                    break;
                case 4:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2) |
                                      (unpacked[i + 2] << 4) |
                                      (unpacked[i + 3] << 6);
                    break;
                case 3:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2) |
                                      (unpacked[i + 2] << 4);
                    break;
                case 2:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 2);
                    break;
                case 1:
                    packed[outByte] = (unpacked[i]);
                    break;
            }
            return packed;
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED pack3(const SEQUENCE_UNPACKED& unpacked,
                                     const ALPHABET& alphabet) {
            lensq outByte = 0;
            lensq i = 0;
            SEQUENCE_PACKED packed = reserveSpaceForPacked<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
            for (; i + 8 <= unpacked.size(); i += 8) {
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                      (unpacked[i + 3] << 1) |
                                      (unpacked[i + 4] << 4) |
                                      (unpacked[i + 5] << 7);
                packed[outByte + 2] = (unpacked[i + 5] >> 1) |
                                      (unpacked[i + 6] << 2) |
                                      (unpacked[i + 7] << 5);
                outByte += 3;
            }
            switch (unpacked.size() - i) {
                case 7:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3) |
                                      (unpacked[i + 2] << 6);
                    packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                          (unpacked[i + 3] << 1) |
                                          (unpacked[i + 4] << 4) |
                                          (unpacked[i + 5] << 7);
                    packed[outByte + 2] = (unpacked[i + 5] >> 1) |
                                          (unpacked[i + 6] << 2);
                    break;
                case 6:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3) |
                                      (unpacked[i + 2] << 6);
                    packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                          (unpacked[i + 3] << 1) |
                                          (unpacked[i + 4] << 4) |
                                          (unpacked[i + 5] << 7);
                    packed[outByte + 2] = (unpacked[i + 5] >> 1);
                    break;
                case 5:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3) |
                                      (unpacked[i + 2] << 6);
                    packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                          (unpacked[i + 3] << 1) |
                                          (unpacked[i + 4] << 4);
                    break;
                case 4:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3) |
                                      (unpacked[i + 2] << 6);
                    packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                          (unpacked[i + 3] << 1);
                    break;
                case 3:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3) |
                                      (unpacked[i + 2] << 6);
                    packed[outByte + 1] = (unpacked[i + 2] >> 2);
                    break;
                case 2:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 3);
                    break;
                case 1:
                    packed[outByte] = (unpacked[i]);
                    break;
            }
            return packed;
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED pack4(const SEQUENCE_UNPACKED& unpacked,
                                     const ALPHABET& alphabet) {
            lensq outByte = 0;
            lensq i = 0;
            SEQUENCE_PACKED packed = reserveSpaceForPacked<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
            for (; i + 8 <= unpacked.size(); i += 8) {
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]) |
                                      (unpacked[i + 3] << 4);
                packed[outByte + 2] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 4);
                packed[outByte + 3] = (unpacked[i + 6]) |
                                      (unpacked[i + 7] << 4);
                outByte += 4;
            }
            switch (unpacked.size() - i) {
                case 7:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    packed[outByte + 1] = (unpacked[i + 2]) |
                                          (unpacked[i + 3] << 4);
                    packed[outByte + 2] = (unpacked[i + 4]) |
                                          (unpacked[i + 5] << 4);
                    packed[outByte + 3] = (unpacked[i + 6]);
                    break;
                case 6:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    packed[outByte + 1] = (unpacked[i + 2]) |
                                          (unpacked[i + 3] << 4);
                    packed[outByte + 2] = (unpacked[i + 4]) |
                                          (unpacked[i + 5] << 4);
                    break;
                case 5:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    packed[outByte + 1] = (unpacked[i + 2]) |
                                          (unpacked[i + 3] << 4);
                    packed[outByte + 2] = (unpacked[i + 4]);
                    break;
                case 4:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    packed[outByte + 1] = (unpacked[i + 2]) |
                                          (unpacked[i + 3] << 4);
                    break;
                case 3:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    packed[outByte + 1] = (unpacked[i + 2]);
                    break;
                case 2:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 4);
                    break;
                case 1:
                    packed[outByte] = (unpacked[i]);
                    break;
            }
            return packed;
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED pack5(const SEQUENCE_UNPACKED& unpacked,
                                     const ALPHABET& alphabet) {
            lensq outByte = 0;
            lensq i = 0;
            SEQUENCE_PACKED packed = reserveSpaceForPacked<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
            for (; i + 8 <= unpacked.size(); i += 8) {
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2) |
                                      (unpacked[i + 3] << 7);
                packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                      (unpacked[i + 4] << 4);
                packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                      (unpacked[i + 5] << 1) |
                                      (unpacked[i + 6] << 6);
                packed[outByte + 4] = (unpacked[i + 6] >> 2) |
                                      (unpacked[i + 7] << 3);
                outByte += 5;
            }
            switch (unpacked.size() - i) {
                case 7:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                          (unpacked[i + 2] << 2) |
                                          (unpacked[i + 3] << 7);
                    packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                          (unpacked[i + 4] << 4);
                    packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                          (unpacked[i + 5] << 1) |
                                          (unpacked[i + 6] << 6);
                    packed[outByte + 4] = (unpacked[i + 6] >> 2);
                    break;
                case 6:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                          (unpacked[i + 2] << 2) |
                                          (unpacked[i + 3] << 7);
                    packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                          (unpacked[i + 4] << 4);
                    packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                          (unpacked[i + 5] << 1);
                    break;
                case 5:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                          (unpacked[i + 2] << 2) |
                                          (unpacked[i + 3] << 7);
                    packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                          (unpacked[i + 4] << 4);
                    packed[outByte + 3] = (unpacked[i + 4] >> 4);
                    break;
                case 4:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                          (unpacked[i + 2] << 2) |
                                          (unpacked[i + 3] << 7);
                    packed[outByte + 2] = (unpacked[i + 3] >> 1);
                    break;
                case 3:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                          (unpacked[i + 2] << 2);
                    break;
                case 2:
                    packed[outByte] = (unpacked[i]) |
                                      (unpacked[i + 1] << 5);
                    packed[outByte + 1] = (unpacked[i + 1] >> 3);
                    break;
                case 1:
                    packed[outByte] = (unpacked[i]);
                    break;
            }
            return packed;
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        SEQUENCE_PACKED pack(const SEQUENCE_UNPACKED& unpacked, const ALPHABET& alphabet) {
            switch (getAlphabetSize(alphabet)) {
                case 2:
                    return pack2<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
                case 3:
                    return pack3<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
                case 4:
                    return pack4<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
                case 5:
                    return pack5<SEQUENCE_UNPACKED, SEQUENCE_PACKED, ALPHABET>(unpacked, alphabet);
                default:
                    throw std::invalid_argument("unpacked");
            }
        }
    }
}

#endif //TIDYSQ_INTERNAL_PACK_H
