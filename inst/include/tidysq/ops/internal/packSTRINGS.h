#ifndef TIDYSQ_PACKSTRINGS_H
#define TIDYSQ_PACKSTRINGS_H

#include <stdexcept>
#include "../../types/general.h"
#include "util.h"

namespace tidysq::internal {

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS2(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, STRINGS, INTERNAL_OUT>(unpacked, alphabet);
        letvalue NAValue = getNAValue(alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                              (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                              (match(unpacked[i + 2], alphabet, NAValue) << 4) |
                              (match(unpacked[i + 3], alphabet, NAValue) << 6);
            packed[outByte + 1] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                  (match(unpacked[i + 5], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 6], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 7], alphabet, NAValue) << 6);
            outByte += 2;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 2) |
                                      (match(unpacked[i + 6], alphabet, NAValue) << 4);
                break;
            case 6:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 2);
                break;
            case 5:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 4], alphabet, NAValue));
                break;
            case 4:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 6);
                break;
            case 3:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 4);
                break;
            case 2:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 2);
                break;
            case 1:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS3(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, STRINGS, INTERNAL_OUT>(unpacked, alphabet);
        letvalue NAValue = getNAValue(alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                              (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                              (match(unpacked[i + 2], alphabet, NAValue) << 6);
            packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 1) |
                                  (match(unpacked[i + 4], alphabet, NAValue) << 4) |
                                  (match(unpacked[i + 5], alphabet, NAValue) << 7);
            packed[outByte + 2] = (match(unpacked[i + 5], alphabet, NAValue) >> 1) |
                                  (match(unpacked[i + 6], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 7], alphabet, NAValue) << 5);
            outByte += 3;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 5], alphabet, NAValue) >> 1) |
                                      (match(unpacked[i + 6], alphabet, NAValue) << 2);
                break;
            case 6:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 5], alphabet, NAValue) >> 1);
                break;
            case 5:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4);
                break;
            case 4:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 1);
                break;
            case 3:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 6);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue) >> 2);
                break;
            case 2:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 3);
                break;
            case 1:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS4(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                     const Alphabet<INTERNAL_IN> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, STRINGS, INTERNAL_OUT>(unpacked, alphabet);
        letvalue NAValue = getNAValue(alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                              (match(unpacked[i + 1], alphabet, NAValue) << 4);
            packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue)) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 4);
            packed[outByte + 2] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                  (match(unpacked[i + 5], alphabet, NAValue) << 4);
            packed[outByte + 3] = (match(unpacked[i + 6], alphabet, NAValue)) |
                                  (match(unpacked[i + 7], alphabet, NAValue) << 4);
            outByte += 4;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue)) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 4);
                packed[outByte + 2] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 4);
                packed[outByte + 3] = (match(unpacked[i + 6], alphabet, NAValue));
                break;
            case 6:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue)) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 4);
                packed[outByte + 2] = (match(unpacked[i + 4], alphabet, NAValue)) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 4);
                break;
            case 5:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue)) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 4);
                packed[outByte + 2] = (match(unpacked[i + 4], alphabet, NAValue));
                break;
            case 4:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue)) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 4);
                break;
            case 3:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                packed[outByte + 1] = (match(unpacked[i + 2], alphabet, NAValue));
                break;
            case 2:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 4);
                break;
            case 1:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS5(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                     const Alphabet<INTERNAL_IN> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, STRINGS, INTERNAL_OUT>(unpacked, alphabet);
        letvalue NAValue = getNAValue(alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                              (match(unpacked[i + 1], alphabet, NAValue) << 5);
            packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                  (match(unpacked[i + 2], alphabet, NAValue) << 2) |
                                  (match(unpacked[i + 3], alphabet, NAValue) << 7);
            packed[outByte + 2] = (match(unpacked[i + 3], alphabet, NAValue) >> 1) |
                                  (match(unpacked[i + 4], alphabet, NAValue) << 4);
            packed[outByte + 3] = (match(unpacked[i + 4], alphabet, NAValue) >> 4) |
                                  (match(unpacked[i + 5], alphabet, NAValue) << 1) |
                                  (match(unpacked[i + 6], alphabet, NAValue) << 6);
            packed[outByte + 4] = (match(unpacked[i + 6], alphabet, NAValue) >> 2) |
                                  (match(unpacked[i + 7], alphabet, NAValue) << 3);
            outByte += 5;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                      (match(unpacked[i + 2], alphabet, NAValue) << 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 3], alphabet, NAValue) >> 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4);
                packed[outByte + 3] = (match(unpacked[i + 4], alphabet, NAValue) >> 4) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 1) |
                                      (match(unpacked[i + 6], alphabet, NAValue) << 6);
                packed[outByte + 4] = (match(unpacked[i + 6], alphabet, NAValue) >> 2);
                break;
            case 6:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                      (match(unpacked[i + 2], alphabet, NAValue) << 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 3], alphabet, NAValue) >> 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4);
                packed[outByte + 3] = (match(unpacked[i + 4], alphabet, NAValue) >> 4) |
                                      (match(unpacked[i + 5], alphabet, NAValue) << 1);
                break;
            case 5:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                      (match(unpacked[i + 2], alphabet, NAValue) << 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 3], alphabet, NAValue) >> 1) |
                                      (match(unpacked[i + 4], alphabet, NAValue) << 4);
                packed[outByte + 3] = (match(unpacked[i + 4], alphabet, NAValue) >> 4);
                break;
            case 4:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                      (match(unpacked[i + 2], alphabet, NAValue) << 2) |
                                      (match(unpacked[i + 3], alphabet, NAValue) << 7);
                packed[outByte + 2] = (match(unpacked[i + 3], alphabet, NAValue) >> 1);
                break;
            case 3:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3) |
                                      (match(unpacked[i + 2], alphabet, NAValue) << 2);
                break;
            case 2:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue)) |
                                  (match(unpacked[i + 1], alphabet, NAValue) << 5);
                packed[outByte + 1] = (match(unpacked[i + 1], alphabet, NAValue) >> 3);
                break;
            case 1:
                packed[outByte] = (match(unpacked[i], alphabet, NAValue));
                break;
        }
        return packed;
    }
    

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                       const Alphabet<INTERNAL_IN> &alphabet) {
        switch (getAlphabetSize(alphabet)) {
            case 2:
                return packSTRINGS2<INTERNAL_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 3:
                return packSTRINGS3<INTERNAL_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 4:
                return packSTRINGS4<INTERNAL_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 5:
                return packSTRINGS5<INTERNAL_IN, INTERNAL_OUT>(unpacked, alphabet);
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACKSTRINGS_H
