#ifndef TIDYSQ_PACKSTRING_STANDARD_H
#define TIDYSQ_PACKSTRING_STANDARD_H

#include <stdexcept>

#include "../../types/general.h"
#include "../../types/SequenceProtoRCPP.h"
#include "../../types/AlphabetRCPP.h"
#include "util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD2(const SequenceProto<STD, STRING> &unpacked,
                                                const Alphabet<STD> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<STD, STRING, INTERNAL_OUT>(unpacked, alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 6);
            packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 7], alphabet) << 6);
            outByte += 2;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 4);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 2);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet));
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 6);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 4);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 2);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD3(const SequenceProto<STD, STRING> &unpacked,
                                       const Alphabet<STD> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<STD, STRING, INTERNAL_OUT>(unpacked, alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
            packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 1) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 7);
            packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) >> 1) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 7], alphabet) << 5);
            outByte += 3;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) >> 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 2);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) >> 1);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4);
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 1);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 6);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) >> 2);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 3);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD4(const SequenceProto<STD, STRING> &unpacked,
                                       const Alphabet<STD> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<STD, STRING, INTERNAL_OUT>(unpacked, alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
            packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 4);
            packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 4);
            packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 7], alphabet) << 4);
            outByte += 4;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet));
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 4);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 4);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet));
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet)) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 4);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet));
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 4);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet));
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD5(const SequenceProto<STD, STRING> &unpacked,
                                       const Alphabet<STD> &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<STD, STRING, INTERNAL_OUT>(unpacked, alphabet);

        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                              (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
            packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 7);
            packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) >> 1) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4);
            packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) >> 4) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 1) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 6);
            packed[outByte + 4] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) >> 2) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 7], alphabet) << 3);
            outByte += 5;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) >> 4) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) << 6);
                packed[outByte + 4] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 6], alphabet) >> 2);
                break;
            case 6:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) >> 4) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 5], alphabet) << 1);
                break;
            case 5:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) >> 1) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) << 4);
                packed[outByte + 3] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 4], alphabet) >> 4);
                break;
            case 4:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) << 7);
                packed[outByte + 2] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 3], alphabet) >> 1);
                break;
            case 3:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3) |
                                      (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 2], alphabet) << 2);
                break;
            case 2:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet)) |
                                  (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) << 5);
                packed[outByte + 1] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i + 1], alphabet) >> 3);
                break;
            case 1:
                packed[outByte] = (ValueToLetterMatcher<STD>::matchStandard(unpacked[i], alphabet));
                break;
        }
        return packed;
    }


    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD(const SequenceProto<STD, STRING> &unpacked,
                                      const Alphabet<STD> &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 2:
                return packSTRING_STANDARD2<INTERNAL_OUT>(unpacked, alphabet);
            case 3:
                return packSTRING_STANDARD3<INTERNAL_OUT>(unpacked, alphabet);
            case 4:
                return packSTRING_STANDARD4<INTERNAL_OUT>(unpacked, alphabet);
            case 5:
                return packSTRING_STANDARD5<INTERNAL_OUT>(unpacked, alphabet);
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }

    template<InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING_STANDARD(const SequenceProto<RCPP, STRING> &unpacked,
                                               const Alphabet<RCPP> &alphabet) {

        return packSTRING_STANDARD<INTERNAL_OUT>(static_cast<SequenceProto<STD, STRING>>(unpacked),
                                                 static_cast<Alphabet<STD>>(alphabet));
    }
}

#endif //TIDYSQ_PACKSTRING_STANDARD_H
