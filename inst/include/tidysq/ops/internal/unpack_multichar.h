#ifndef TIDYSQ_UNPACK_MULTICHAR_H
#define TIDYSQ_UNPACK_MULTICHAR_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    void unpack_multichar_string_2(const Sequence<INTERNAL_IN> &packed,
                                   ProtoSequence<INTERNAL_OUT, STRING> &unpacked,
                                   const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = packed.originalLength();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte] >> 6) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1]) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 3, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 6) & 3, alphabet);
            in_byte += 2;
        }
        switch (out_len - i) {
            case 7:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 6) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 3, alphabet);
                break;
            case 6:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 6) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 3, alphabet);
                break;
            case 5:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 6) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 3, alphabet);
                break;
            case 4:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 6) & 3, alphabet);
                break;
            case 3:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 3, alphabet);
                break;
            case 2:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 2) & 3, alphabet);
                break;
            case 1:
                unpacked += matchLetterMultichar((packed[in_byte]) & 3, alphabet);
                break;
        }
    }


    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    void unpack_multichar_string_3(const Sequence<INTERNAL_IN> &packed,
                                   ProtoSequence<INTERNAL_OUT, STRING> &unpacked,
                                   const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = packed.originalLength();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                     ((packed[in_byte + 1] << 2) & 7), alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 1) & 7, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 7, alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                     ((packed[in_byte + 2] << 1) & 7), alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 2] >> 2) & 7, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 2] >> 5) & 7, alphabet);
            in_byte += 3;
        }
        switch (out_len - i) {
            case 7:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 7), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2] >> 2) & 7, alphabet);
                break;
            case 6:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 7), alphabet);
                break;
            case 5:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 7, alphabet);
                break;
            case 4:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 1) & 7, alphabet);
                break;
            case 3:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                break;
            case 2:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 3) & 7, alphabet);
                break;
            case 1:
                unpacked += matchLetterMultichar((packed[in_byte]) & 7, alphabet);
                break;
        }
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    void unpack_multichar_string_4(const Sequence<INTERNAL_IN> &packed,
                                   ProtoSequence<INTERNAL_OUT, STRING> &unpacked,
                                   const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = packed.originalLength();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 2]) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 2] >> 4) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 3]) & 15, alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 3] >> 4) & 15, alphabet);
            in_byte += 4;
        }
        switch (out_len - i) {
            case 7:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 3]) & 15, alphabet);
                break;
            case 6:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2] >> 4) & 15, alphabet);
                break;
            case 5:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 2]) & 15, alphabet);
                break;
            case 4:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 4) & 15, alphabet);
                break;
            case 3:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1]) & 15, alphabet);
                break;
            case 2:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                unpacked += matchLetterMultichar((packed[in_byte] >> 4) & 15, alphabet);
                break;
            case 1:
                unpacked += matchLetterMultichar((packed[in_byte]) & 15, alphabet);
                break;
        }
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    void unpack_multichar_string_5(const Sequence<INTERNAL_IN> &packed,
                                   ProtoSequence<INTERNAL_OUT, STRING> &unpacked,
                                   const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = packed.originalLength();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                     ((packed[in_byte + 1] << 3) & 31), alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                     ((packed[in_byte + 2] << 1) & 31), alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte + 2] >> 4) & 15) |
                                                     ((packed[in_byte + 3] << 4) & 31), alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 3] >> 1) & 31, alphabet);
            unpacked += matchLetterMultichar(((packed[in_byte + 3] >> 6) & 3) |
                                                     ((packed[in_byte + 4] << 2) & 31), alphabet);
            unpacked += matchLetterMultichar((packed[in_byte + 4] >> 3) & 31, alphabet);
            in_byte += 5;
        }
        switch (out_len - i) {
            case 7:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 2] >> 4) & 15) |
                                                         ((packed[in_byte + 3] << 4) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 3] >> 1) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 3] >> 6) & 3) |
                                                         ((packed[in_byte + 4] << 2) & 31), alphabet);
                break;
            case 6:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 2] >> 4) & 15) |
                                                         ((packed[in_byte + 3] << 4) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 3] >> 1) & 31, alphabet);
                break;
            case 5:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 2] >> 4) & 15) |
                                                         ((packed[in_byte + 3] << 4) & 31), alphabet);
                break;
            case 4:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 31), alphabet);
                break;
            case 3:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked += matchLetterMultichar((packed[in_byte + 1] >> 2) & 31, alphabet);
                break;
            case 2:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                unpacked += matchLetterMultichar(((packed[in_byte] >> 5) & 7) |
                                                         ((packed[in_byte + 1] << 3) & 31), alphabet);
                break;
            case 1:
                unpacked += matchLetterMultichar((packed[in_byte]) & 31, alphabet);
                break;
        }
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    void unpack_multichar_string(const Sequence<INTERNAL_IN> &packed,
                                 ProtoSequence<INTERNAL_OUT, STRING> &unpacked,
                                 const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                unpack_multichar_string_2(packed, unpacked, alphabet);
                return;
            case 3:
                unpack_multichar_string_3(packed, unpacked, alphabet);
                return;
            case 4:
                unpack_multichar_string_4(packed, unpacked, alphabet);
                return;
            case 5:
                unpack_multichar_string_5(packed, unpacked, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_UNPACK_MULTICHAR_H
