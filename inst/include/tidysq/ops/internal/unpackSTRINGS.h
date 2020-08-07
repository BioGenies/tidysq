#ifndef TIDYSQ_UNPACKSTRINGS_H
#define TIDYSQ_UNPACKSTRINGS_H

#include <stdexcept>

#include "../../types/general.h"
#include "util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    SequenceProto<INTERNAL_OUT, STRINGS> unpackSTRINGS2(const Sequence<INTERNAL_IN> &packed,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<INTERNAL_OUT, STRINGS> unpacked = reserveSpaceForUnpacked<INTERNAL_IN, INTERNAL_OUT, STRINGS>(
                packed);
        lensq out_len = unpacked.size();

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
            unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
            unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
            unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 6) & 3, alphabet);
            unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 3, alphabet);
            unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 3, alphabet);
            unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 3, alphabet);
            unpacked[i + 7] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 6) & 3, alphabet);
            in_byte += 2;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 6) & 3, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 3, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 3, alphabet);
                unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 3, alphabet);
                break;
            case 6:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 6) & 3, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 3, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 3, alphabet);
                break;
            case 5:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 6) & 3, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 3, alphabet);
                break;
            case 4:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 6) & 3, alphabet);
                break;
            case 3:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 3, alphabet);
                break;
            case 2:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 2) & 3, alphabet);
                break;
            case 1:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 3, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    SequenceProto<INTERNAL_OUT, STRINGS> unpackSTRINGS3(const Sequence<INTERNAL_IN> &packed,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<INTERNAL_OUT, STRINGS> unpacked = reserveSpaceForUnpacked<INTERNAL_IN, INTERNAL_OUT, STRINGS>(
                packed);
        lensq out_len = unpacked.size();

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
            unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
            unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                       ((packed[in_byte + 1] << 2) & 7), alphabet);
            unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 1) & 7, alphabet);
            unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 7, alphabet);
            unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                       ((packed[in_byte + 2] << 1) & 7), alphabet);
            unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 2) & 7, alphabet);
            unpacked[i + 7] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 5) & 7, alphabet);
            in_byte += 3;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                           ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 7), alphabet);
                unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 2) & 7, alphabet);
                break;
            case 6:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                           ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 7, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 7), alphabet);
                break;
            case 5:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                           ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 1) & 7, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 7, alphabet);
                break;
            case 4:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                           ((packed[in_byte + 1] << 2) & 7), alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 1) & 7, alphabet);
                break;
            case 3:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 6) & 3) |
                                           ((packed[in_byte + 1] << 2) & 7), alphabet);
                break;
            case 2:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 3) & 7, alphabet);
                break;
            case 1:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 7, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    SequenceProto<INTERNAL_OUT, STRINGS> unpackSTRINGS4(const Sequence<INTERNAL_IN> &packed,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<INTERNAL_OUT, STRINGS> unpacked = reserveSpaceForUnpacked<INTERNAL_IN, INTERNAL_OUT, STRINGS>(
                packed);
        lensq out_len = unpacked.size();

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
            unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
            unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
            unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 15, alphabet);
            unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2]) & 15, alphabet);
            unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 4) & 15, alphabet);
            unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3]) & 15, alphabet);
            unpacked[i + 7] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3] >> 4) & 15, alphabet);
            in_byte += 4;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2]) & 15, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 4) & 15, alphabet);
                unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3]) & 15, alphabet);
                break;
            case 6:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2]) & 15, alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2] >> 4) & 15, alphabet);
                break;
            case 5:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 15, alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 2]) & 15, alphabet);
                break;
            case 4:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 4) & 15, alphabet);
                break;
            case 3:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1]) & 15, alphabet);
                break;
            case 2:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte] >> 4) & 15, alphabet);
                break;
            case 1:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 15, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    SequenceProto<INTERNAL_OUT, STRINGS> unpackSTRINGS5(const Sequence<INTERNAL_IN> &packed,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        lensq in_byte = 0;
        SequenceProto<INTERNAL_OUT, STRINGS> unpacked = reserveSpaceForUnpacked<INTERNAL_IN, INTERNAL_OUT, STRINGS>(
                packed);
        lensq out_len = unpacked.size();

        lensq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
            unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                       ((packed[in_byte + 1] << 3) & 31), alphabet);
            unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
            unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                       ((packed[in_byte + 2] << 1) & 31), alphabet);
            unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 2] >> 4) & 15) |
                                       ((packed[in_byte + 3] << 4) & 31), alphabet);
            unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3] >> 1) & 31, alphabet);
            unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 3] >> 6) & 3) |
                                       ((packed[in_byte + 4] << 2) & 31), alphabet);
            unpacked[i + 7] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 4] >> 3) & 31, alphabet);
            in_byte += 5;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 2] >> 4) & 15) |
                                           ((packed[in_byte + 3] << 4) & 31), alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3] >> 1) & 31, alphabet);
                unpacked[i + 6] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 3] >> 6) & 3) |
                                           ((packed[in_byte + 4] << 2) & 31), alphabet);
                break;
            case 6:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 2] >> 4) & 15) |
                                           ((packed[in_byte + 3] << 4) & 31), alphabet);
                unpacked[i + 5] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 3] >> 1) & 31, alphabet);
                break;
            case 5:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 31), alphabet);
                unpacked[i + 4] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 2] >> 4) & 15) |
                                           ((packed[in_byte + 3] << 4) & 31), alphabet);
                break;
            case 4:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
                unpacked[i + 3] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte + 1] >> 7) & 1) |
                                           ((packed[in_byte + 2] << 1) & 31), alphabet);
                break;
            case 3:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                unpacked[i + 2] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte + 1] >> 2) & 31, alphabet);
                break;
            case 2:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                unpacked[i + 1] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match(((packed[in_byte] >> 5) & 7) |
                                           ((packed[in_byte + 1] << 3) & 31), alphabet);
                break;
            case 1:
                unpacked[i] = LetterToValueMatcher<INTERNAL_IN, INTERNAL_OUT>::match((packed[in_byte]) & 31, alphabet);
                break;
        }
        return unpacked;
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    SequenceProto<INTERNAL_OUT, STRINGS> unpackSTRINGS(const Sequence<INTERNAL_IN> &packed,
                                                       const Alphabet<INTERNAL_IN> &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 2:
                return unpackSTRINGS2<INTERNAL_IN, INTERNAL_OUT>(packed, alphabet);
            case 3:
                return unpackSTRINGS3<INTERNAL_IN, INTERNAL_OUT>(packed, alphabet);
            case 4:
                return unpackSTRINGS4<INTERNAL_IN, INTERNAL_OUT>(packed, alphabet);
            case 5:
                return unpackSTRINGS5<INTERNAL_IN, INTERNAL_OUT>(packed, alphabet);
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_UNPACKSTRINGS_H
