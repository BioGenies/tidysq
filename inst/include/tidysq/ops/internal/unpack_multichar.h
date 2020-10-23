#ifndef TIDYSQ_UNPACK_MULTICHAR_H
#define TIDYSQ_UNPACK_MULTICHAR_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    void unpack_multichar_2(const Sequence<INTERNAL_IN> &packed,
                            ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                            const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.length();
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

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    void unpack_multichar(const Sequence<INTERNAL_IN> &packed,
                          ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                          const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                unpack_multichar_2(packed, unpacked, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_UNPACK_MULTICHAR_H
