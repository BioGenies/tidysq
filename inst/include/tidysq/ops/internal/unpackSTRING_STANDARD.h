#ifndef TIDYSQ_UNPACKSTRING_STANDARD_H
#define TIDYSQ_UNPACKSTRING_STANDARD_H

#include <stdexcept>

#include "../../types/general.h"
#include "util.h"
#include "../../types/SequenceProtoSTD.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STANDARD2(const Sequence <INTERNAL_IN> &packed,
                                                      const Alphabet<STD> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<STD, STRING> unpacked = SequenceProto<STD, STRING>();
        lensq out_len = getOriginalLength(packed);
        unpacked.reserve(out_len);

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += {
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 6) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 6) & 3, alphabet)};
            in_byte += 2;
        }
        switch (out_len - i) {
            case 7:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 6) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 3, alphabet)
                };
                break;
            case 6:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 6) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 3, alphabet)
                };
                break;
            case 5:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 6) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 3, alphabet)
                };
                break;
            case 4:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 6) & 3, alphabet)
                };
                break;
            case 3:
                unpacked += {
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 3, alphabet)
                };
                break;
            case 2:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 2) & 3, alphabet)
                };
                break;
            case 1:
                unpacked += LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 3, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STANDARD3(const Sequence<INTERNAL_IN> &packed,
                                                      const Alphabet<STD> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<STD, STRING> unpacked = SequenceProto<STD, STRING>();
        lensq out_len = getOriginalLength(packed);
        unpacked.reserve(out_len);

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += {
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                    ((packed[in_byte + 1] << 2) & 7), alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 1) & 7, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 7, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                    ((packed[in_byte + 2] << 1) & 7), alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 2) & 7, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 5) & 7, alphabet)
            };
            in_byte += 3;
        }
        switch (out_len - i) {
            case 7:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                        ((packed[in_byte + 1] << 2) & 7),
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 1) & 7,
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 7,
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                        ((packed[in_byte + 2] << 1) & 7),
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 2) & 7,
                        alphabet)
                };
                break;
            case 6:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                        ((packed[in_byte + 1] << 2) & 7),
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 1) & 7,
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 7,
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                        ((packed[in_byte + 2] << 1) & 7),
                        alphabet)
                };
                break;
            case 5:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                        ((packed[in_byte + 1] << 2) & 7),
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 1) & 7,
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 7,
                        alphabet)
                };
                break;
            case 4:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                        ((packed[in_byte + 1] << 2) & 7),
                        alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 1) & 7,
                        alphabet)
                };
                break;
            case 3:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 6) & 3) |
                        ((packed[in_byte + 1] << 2) & 7),
                        alphabet)
                };
                break;
            case 2:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 3) & 7, alphabet)
                };
                break;
            case 1:
                unpacked += LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 7, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STANDARD4(const Sequence <INTERNAL_IN> &packed,
                                                         const Alphabet<STD> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<STD, STRING> unpacked = SequenceProto<STD, STRING>();
        lensq out_len = getOriginalLength(packed);
        unpacked.reserve(out_len);

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += {
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2]) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 4) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3]) & 15, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3] >> 4) & 15, alphabet)
            };
            in_byte += 4;
        }
        switch (out_len - i) {
            case 7:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 15,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 4) & 15,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3]) & 15, alphabet)
                };
                break;
            case 6:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 15,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2] >> 4) & 15,
                                                                      alphabet)
                };
                break;
            case 5:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 15,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 2]) & 15, alphabet)
                };
                break;
            case 4:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 4) & 15,
                                                                      alphabet)
                };
                break;
            case 3:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1]) & 15, alphabet)
                };
                break;
            case 2:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte] >> 4) & 15, alphabet)
                };
                break;
            case 1:
                unpacked += LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 15, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STANDARD5(const Sequence <INTERNAL_IN> &packed,
                                                         const Alphabet<STD> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<STD, STRING> unpacked = SequenceProto<STD, STRING>();
        lensq out_len = getOriginalLength(packed);
        unpacked.reserve(out_len);

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked += {
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                  ((packed[in_byte + 1] << 3) & 31),
                                                                  alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                                                                  ((packed[in_byte + 2] << 1) & 31),
                                                                  alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 2] >> 4) & 15) |
                                                                  ((packed[in_byte + 3] << 4) & 31),
                                                                  alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3] >> 1) & 31, alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 3] >> 6) & 3) |
                                                                  ((packed[in_byte + 4] << 2) & 31),
                                                                  alphabet),
                    LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 4] >> 3) & 31, alphabet)
            };
            in_byte += 5;
        }
        switch (out_len - i) {
            case 7:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                                                                      ((packed[in_byte + 2] << 1) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 2] >> 4) & 15) |
                                                                      ((packed[in_byte + 3] << 4) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3] >> 1) & 31,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 3] >> 6) & 3) |
                                                                      ((packed[in_byte + 4] << 2) & 31),
                                                                      alphabet)
                };
                break;
            case 6:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                                                                      ((packed[in_byte + 2] << 1) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 2] >> 4) & 15) |
                                                                      ((packed[in_byte + 3] << 4) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 3] >> 1) & 31,
                                                                      alphabet)
                };
                break;
            case 5:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                                                                      ((packed[in_byte + 2] << 1) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 2] >> 4) & 15) |
                                                                      ((packed[in_byte + 3] << 4) & 31),
                                                                      alphabet)
                };
                break;
            case 4:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31,
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte + 1] >> 7) & 1) |
                                                                      ((packed[in_byte + 2] << 1) & 31),
                                                                      alphabet)
                };
                break;
            case 3:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte + 1] >> 2) & 31,
                                                                      alphabet)
                };
                break;
            case 2:
                unpacked += {
                        LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet),
                        LetterToValueMatcher<STD, STD>::matchStandard(((packed[in_byte] >> 5) & 7) |
                                                                      ((packed[in_byte + 1] << 3) & 31),
                                                                      alphabet)
                };
                break;
            case 1:
                unpacked += LetterToValueMatcher<STD, STD>::matchStandard((packed[in_byte]) & 31, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STANDARD(const Sequence <INTERNAL_IN> &packed,
                                                       const Alphabet<STD> &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 2:
                return unpackSTRING_STANDARD2<INTERNAL_IN>(packed, alphabet);
            case 3:
                return unpackSTRING_STANDARD3<INTERNAL_IN>(packed, alphabet);
            case 4:
                return unpackSTRING_STANDARD4<INTERNAL_IN>(packed, alphabet);
            case 5:
                return unpackSTRING_STANDARD5<INTERNAL_IN>(packed, alphabet);
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_UNPACKSTRING_STANDARD_H
