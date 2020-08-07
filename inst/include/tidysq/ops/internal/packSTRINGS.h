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

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 6);
            packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 7], alphabet) << 6);
            outByte += 2;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 4);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 2);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet));
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 6);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 4);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 2);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet));
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

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
            packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 1) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 7);
            packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) >> 1) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 7], alphabet) << 5);
            outByte += 3;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) >> 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 2);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) >> 1);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4);
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 1);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) >> 2);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 3);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet));
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

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
            packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 4);
            packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 4);
            packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 7], alphabet) << 4);
            outByte += 4;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet));
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 4);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet));
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 4);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet));
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 4);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet));
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

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
            packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 7);
            packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) >> 1) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4);
            packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) >> 4) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 1) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 6);
            packed[outByte + 4] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) >> 2) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 7], alphabet) << 3);
            outByte += 5;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) >> 4) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) << 6);
                packed[outByte + 4] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 6], alphabet) >> 2);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) >> 4) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 5], alphabet) << 1);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 4], alphabet) >> 4);
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 3], alphabet) >> 1);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 2], alphabet) << 2);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i + 1], alphabet) >> 3);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<INTERNAL_IN>::match(unpacked[i], alphabet));
                break;
        }
        return packed;
    }
    

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRINGS(const SequenceProto<INTERNAL_IN, STRINGS> &unpacked,
                                       const Alphabet<INTERNAL_IN> &alphabet) {
        switch (alphabet.alphabetSize()) {
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
