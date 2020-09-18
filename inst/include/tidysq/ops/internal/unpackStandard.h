#ifndef TIDYSQ_UNPACKSTANDARD_H
#define TIDYSQ_UNPACKSTANDARD_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/SequenceProto.h"
#include "tidysq/ops/internal/utilp.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    void unpackStandard3(const Sequence<INTERNAL_IN> &packed,
                         SequenceProto<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        lensq in_byte = 0;
        lensq out_len = unpacked.size();
        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
            unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
            unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                     ((packed[in_byte + 1] << 2) & 7), alphabet);
            unpacked[i + 3] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 1) & 7, alphabet);
            unpacked[i + 4] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 4) & 7, alphabet);
            unpacked[i + 5] = matchLetter<PROTO_OUT>(((packed[in_byte + 1] >> 7) & 1) |
                                                     ((packed[in_byte + 2] << 1) & 7), alphabet);
            unpacked[i + 6] = matchLetter<PROTO_OUT>((packed[in_byte + 2] >> 2) & 7, alphabet);
            unpacked[i + 7] = matchLetter<PROTO_OUT>((packed[in_byte + 2] >> 5) & 7, alphabet);
            in_byte += 3;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked[i + 5] = matchLetter<PROTO_OUT>(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 7), alphabet);
                unpacked[i + 6] = matchLetter<PROTO_OUT>((packed[in_byte + 2] >> 2) & 7, alphabet);
                break;
            case 6:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked[i + 5] = matchLetter<PROTO_OUT>(((packed[in_byte + 1] >> 7) & 1) |
                                                         ((packed[in_byte + 2] << 1) & 7), alphabet);
                break;
            case 5:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 4) & 7, alphabet);
                break;
            case 4:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = matchLetter<PROTO_OUT>((packed[in_byte + 1] >> 1) & 7, alphabet);
                break;
            case 3:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = matchLetter<PROTO_OUT>(((packed[in_byte] >> 6) & 3) |
                                                         ((packed[in_byte + 1] << 2) & 7), alphabet);
                break;
            case 2:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = matchLetter<PROTO_OUT>((packed[in_byte] >> 3) & 7, alphabet);
                break;
            case 1:
                unpacked[i] = matchLetter<PROTO_OUT>((packed[in_byte]) & 7, alphabet);
                break;
        }
    }


    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    void unpackStandard(const Sequence<INTERNAL_IN> &packed,
                        SequenceProto<INTERNAL_OUT, PROTO_OUT> &unpacked,
                        const Alphabet &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 3:
                unpackStandard3(packed, unpacked, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_UNPACKSTANDARD_H
