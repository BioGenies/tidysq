#ifndef TIDYSQ_PACKSTANDARD_H
#define TIDYSQ_PACKSTANDARD_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void packStandard3(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq outByte = 0;
        LenSq i = 0;
        for (; i + 8 <= unpacked.length(); i += 8) {
            packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                              (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                              (unpacked.getLetterValue(i + 2, alphabet) << 6u);
            packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 1u) |
                                  (unpacked.getLetterValue(i + 4, alphabet) << 4u) |
                                  (unpacked.getLetterValue(i + 5, alphabet) << 7u);
            packed[outByte + 2] = (unpacked.getLetterValue(i + 5, alphabet) >> 1u) |
                                  (unpacked.getLetterValue(i + 6, alphabet) << 2u) |
                                  (unpacked.getLetterValue(i + 7, alphabet) << 5u);
            outByte += 3;
        }
        switch (unpacked.length() - i) {
            case 7:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 6u);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 1u) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4u) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 7u);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 5, alphabet) >> 1u) |
                                      (unpacked.getLetterValue(i + 6, alphabet) << 2u);
                break;
            case 6:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 6u);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 1u) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4u) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 7u);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 5, alphabet) >> 1u);
                break;
            case 5:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 6u);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 1u) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4u);
                break;
            case 4:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 6u);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 1u);
                break;
            case 3:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 6u);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet) >> 2u);
                break;
            case 2:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 3u);
                break;
            case 1:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet));
                break;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void packStandard(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                      Sequence<INTERNAL_OUT> &packed,
                      const Alphabet &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 3:
                packStandard3(unpacked, packed, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACKSTANDARD_H
