#ifndef TIDYSQ_PACK_SIMPLE_H
#define TIDYSQ_PACK_SIMPLE_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void pack2(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq outByte = 0;
        LenSq i = 0;
        for (; i + 8 <= unpacked.length(); i += 8) {
            packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                              (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                              (unpacked.getLetterValue(i + 2, alphabet) << 4) |
                              (unpacked.getLetterValue(i + 3, alphabet) << 6);
            packed[outByte + 1] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                  (unpacked.getLetterValue(i + 5, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 6, alphabet) << 4) |
                                  (unpacked.getLetterValue(i + 7, alphabet) << 6);
            outByte += 2;
        }
        switch (unpacked.length() - i) {
            case 7:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 4) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 6);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 2) |
                                      (unpacked.getLetterValue(i + 6, alphabet) << 4);
                break;
            case 6:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 4) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 6);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 2);
                break;
            case 5:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 4) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 6);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 4, alphabet));
                break;
            case 4:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 4) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 6);
                break;
            case 3:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 4);
                break;
            case 2:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 2);
                break;
            case 1:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet));
                break;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void pack3(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
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
    void pack4(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq outByte = 0;
        LenSq i = 0;
        for (; i + 8 <= unpacked.length(); i += 8) {
            packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                              (unpacked.getLetterValue(i + 1, alphabet) << 4);
            packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet)) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 4);
            packed[outByte + 2] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                  (unpacked.getLetterValue(i + 5, alphabet) << 4);
            packed[outByte + 3] = (unpacked.getLetterValue(i + 6, alphabet)) |
                                  (unpacked.getLetterValue(i + 7, alphabet) << 4);
            outByte += 4;
        }
        switch (unpacked.length() - i) {
            case 7:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet)) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 4);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 4);
                packed[outByte + 3] = (unpacked.getLetterValue(i + 6, alphabet));
                break;
            case 6:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet)) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 4);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 4, alphabet)) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 4);
                break;
            case 5:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet)) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 4);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 4, alphabet));
                break;
            case 4:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet)) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 4);
                break;
            case 3:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 2, alphabet));
                break;
            case 2:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 4);
                break;
            case 1:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet));
                break;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void pack5(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq outByte = 0;
        LenSq i = 0;
        for (; i + 8 <= unpacked.length(); i += 8) {
            packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                              (unpacked.getLetterValue(i + 1, alphabet) << 5);
            packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                  (unpacked.getLetterValue(i + 2, alphabet) << 2) |
                                  (unpacked.getLetterValue(i + 3, alphabet) << 7);
            packed[outByte + 2] = (unpacked.getLetterValue(i + 3, alphabet) >> 1) |
                                  (unpacked.getLetterValue(i + 4, alphabet) << 4);
            packed[outByte + 3] = (unpacked.getLetterValue(i + 4, alphabet) >> 4) |
                                  (unpacked.getLetterValue(i + 5, alphabet) << 1) |
                                  (unpacked.getLetterValue(i + 6, alphabet) << 6);
            packed[outByte + 4] = (unpacked.getLetterValue(i + 6, alphabet) >> 2) |
                                  (unpacked.getLetterValue(i + 7, alphabet) << 3);
            outByte += 5;
        }
        switch (unpacked.length() - i) {
            case 7:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                      (unpacked.getLetterValue(i + 2, alphabet) << 2) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 7);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 3, alphabet) >> 1) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4);
                packed[outByte + 3] = (unpacked.getLetterValue(i + 4, alphabet) >> 4) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 1) |
                                      (unpacked.getLetterValue(i + 6, alphabet) << 6);
                packed[outByte + 4] = (unpacked.getLetterValue(i + 6, alphabet) >> 2);
                break;
            case 6:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                      (unpacked.getLetterValue(i + 2, alphabet) << 2) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 7);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 3, alphabet) >> 1) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4);
                packed[outByte + 3] = (unpacked.getLetterValue(i + 4, alphabet) >> 4) |
                                      (unpacked.getLetterValue(i + 5, alphabet) << 1);
                break;
            case 5:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                      (unpacked.getLetterValue(i + 2, alphabet) << 2) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 7);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 3, alphabet) >> 1) |
                                      (unpacked.getLetterValue(i + 4, alphabet) << 4);
                packed[outByte + 3] = (unpacked.getLetterValue(i + 4, alphabet) >> 4);
                break;
            case 4:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                      (unpacked.getLetterValue(i + 2, alphabet) << 2) |
                                      (unpacked.getLetterValue(i + 3, alphabet) << 7);
                packed[outByte + 2] = (unpacked.getLetterValue(i + 3, alphabet) >> 1);
                break;
            case 3:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3) |
                                      (unpacked.getLetterValue(i + 2, alphabet) << 2);
                break;
            case 2:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet)) |
                                  (unpacked.getLetterValue(i + 1, alphabet) << 5);
                packed[outByte + 1] = (unpacked.getLetterValue(i + 1, alphabet) >> 3);
                break;
            case 1:
                packed[outByte] = (unpacked.getLetterValue(i, alphabet));
                break;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void pack_simple(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                     Sequence<INTERNAL_OUT> &packed,
                     const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                pack2(unpacked, packed, alphabet);
                return;
            case 3:
                pack3(unpacked, packed, alphabet);
                return;
            case 4:
                pack4(unpacked, packed, alphabet);
                return;
            case 5:
                pack5(unpacked, packed, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACK_SIMPLE_H
