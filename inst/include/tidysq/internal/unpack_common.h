#pragma once

#include <stdexcept>

#include "tidysq/Alphabet.h"
#include "tidysq/Sequence.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/util/calculate_length.h"
#include "tidysq/util/match_letter.h"

namespace tidysq::internal {
    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common_2(const Sequence<INTERNAL_IN> &packed,
                         ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.size();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 3, alphabet);
            unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 2) & 3, alphabet);
            unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 4) & 3, alphabet);
            unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 6) & 3, alphabet);
            unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 1)) & 3, alphabet);
            unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 2) & 3, alphabet);
            unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 3, alphabet);
            unpacked[i + 7] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 6) & 3, alphabet);
            in_byte += 2;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 3, alphabet);
            case 6:
                unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 2) & 3, alphabet);
            case 5:
                unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 1)) & 3, alphabet);
            case 4:
                unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 6) & 3, alphabet);
            case 3:
                unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 4) & 3, alphabet);
            case 2:
                unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 2) & 3, alphabet);
            case 1:
                unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 3, alphabet);
        }
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common_3(const Sequence<INTERNAL_IN> &packed,
                         ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.size();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 7, alphabet);
            unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 3) & 7, alphabet);
            unpacked[i + 2] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 6) & 3) |
                                                            ((packed(in_byte + 1) << 2) & 7), alphabet);
            unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 1) & 7, alphabet);
            unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 7, alphabet);
            unpacked[i + 5] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 7) & 1) |
                                                            ((packed(in_byte + 2) << 1) & 7), alphabet);
            unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 2) & 7, alphabet);
            unpacked[i + 7] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 5) & 7, alphabet);
            in_byte += 3;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 2) & 7, alphabet);
            case 6:
                unpacked[i + 5] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 7) & 1) |
                                                                ((packed(in_byte + 2) << 1) & 7), alphabet);
            case 5:
                unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 7, alphabet);
            case 4:
                unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 1) & 7, alphabet);
            case 3:
                unpacked[i + 2] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 6) & 3) |
                                                                ((packed(in_byte + 1) << 2) & 7), alphabet);
            case 2:
                unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 3) & 7, alphabet);
            case 1:
                unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 7, alphabet);
        }
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common_4(const Sequence<INTERNAL_IN> &packed,
                         ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.size();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 15, alphabet);
            unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 4) & 15, alphabet);
            unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte + 1)) & 15, alphabet);
            unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 15, alphabet);
            unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 2)) & 15, alphabet);
            unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 4) & 15, alphabet);
            unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 3)) & 15, alphabet);
            unpacked[i + 7] = util::match_letter<PROTO_OUT>((packed(in_byte + 3) >> 4) & 15, alphabet);
            in_byte += 4;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i + 6] = util::match_letter<PROTO_OUT>((packed(in_byte + 3)) & 15, alphabet);
            case 6:
                unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 4) & 15, alphabet);
            case 5:
                unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 2)) & 15, alphabet);
            case 4:
                unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 4) & 15, alphabet);
            case 3:
                unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte + 1)) & 15, alphabet);
            case 2:
                unpacked[i + 1] = util::match_letter<PROTO_OUT>((packed(in_byte) >> 4) & 15, alphabet);
            case 1:
                unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 15, alphabet);
        }
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common_5(const Sequence<INTERNAL_IN> &packed,
                         ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.size();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 31, alphabet);
            unpacked[i + 1] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 5) & 7) |
                                                            ((packed(in_byte + 1) << 3) & 31), alphabet);
            unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 2) & 31, alphabet);
            unpacked[i + 3] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 7) & 1) |
                                                            ((packed(in_byte + 2) << 1) & 31), alphabet);
            unpacked[i + 4] = util::match_letter<PROTO_OUT>(((packed(in_byte + 2) >> 4) & 15) |
                                                            ((packed(in_byte + 3) << 4) & 31), alphabet);
            unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 3) >> 1) & 31, alphabet);
            unpacked[i + 6] = util::match_letter<PROTO_OUT>(((packed(in_byte + 3) >> 6) & 3) |
                                                            ((packed(in_byte + 4) << 2) & 31), alphabet);
            unpacked[i + 7] = util::match_letter<PROTO_OUT>((packed(in_byte + 4) >> 3) & 31, alphabet);
            in_byte += 5;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i + 6] = util::match_letter<PROTO_OUT>(((packed(in_byte + 3) >> 6) & 3) |
                                                                ((packed(in_byte + 4) << 2) & 31), alphabet);
            case 6:
                unpacked[i + 5] = util::match_letter<PROTO_OUT>((packed(in_byte + 3) >> 1) & 31, alphabet);
            case 5:
                unpacked[i + 4] = util::match_letter<PROTO_OUT>(((packed(in_byte + 2) >> 4) & 15) |
                                                                ((packed(in_byte + 3) << 4) & 31), alphabet);
            case 4:
                unpacked[i + 3] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 7) & 1) |
                                                                ((packed(in_byte + 2) << 1) & 31), alphabet);
            case 3:
                unpacked[i + 2] = util::match_letter<PROTO_OUT>((packed(in_byte + 1) >> 2) & 31, alphabet);
            case 2:
                unpacked[i + 1] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 5) & 7) |
                                                                ((packed(in_byte + 1) << 3) & 31), alphabet);
            case 1:
                unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 31, alphabet);
        }
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common_6(const Sequence<INTERNAL_IN> &packed,
                         ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                         const Alphabet &alphabet) {
        LenSq in_byte = 0;
        LenSq out_len = unpacked.size();
        LenSq i = 0;
        for (; i + 8 <= out_len; i += 8) {
            unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 63, alphabet);
            unpacked[i + 1] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 6) & 3) |
                                                            ((packed(in_byte + 1) << 2) & 63), alphabet);
            unpacked[i + 2] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 4) & 15) |
                                                            ((packed(in_byte + 2) << 4) & 63), alphabet);
            unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 2) & 63, alphabet);
            unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 3)) & 63, alphabet);
            unpacked[i + 5] = util::match_letter<PROTO_OUT>(((packed(in_byte + 3) >> 6) & 3) |
                                                            ((packed(in_byte + 4) << 2) & 63), alphabet);
            unpacked[i + 6] = util::match_letter<PROTO_OUT>(((packed(in_byte + 4) >> 4) & 15) |
                                                            ((packed(in_byte + 5) << 4) & 63), alphabet);
            unpacked[i + 7] = util::match_letter<PROTO_OUT>((packed(in_byte + 5) >> 2) & 63, alphabet);
            in_byte += 6;
        }
        switch (out_len - i) {
            case 7:
                unpacked[i + 6] = util::match_letter<PROTO_OUT>(((packed(in_byte + 4) >> 4) & 63) |
                                                                ((packed(in_byte + 5) << 4) & 63), alphabet);
            case 6:
                unpacked[i + 5] = util::match_letter<PROTO_OUT>(((packed(in_byte + 3) >> 6) & 15) |
                                                                ((packed(in_byte + 4) << 2) & 63), alphabet);
            case 5:
                unpacked[i + 4] = util::match_letter<PROTO_OUT>((packed(in_byte + 3)) & 63, alphabet);
            case 4:
                unpacked[i + 3] = util::match_letter<PROTO_OUT>((packed(in_byte + 2) >> 2) & 63, alphabet);
            case 3:
                unpacked[i + 2] = util::match_letter<PROTO_OUT>(((packed(in_byte + 1) >> 4) & 15) |
                                                                ((packed(in_byte + 2) << 4) & 63), alphabet);
            case 2:
                unpacked[i + 1] = util::match_letter<PROTO_OUT>(((packed(in_byte) >> 6) & 3) |
                                                                ((packed(in_byte + 1) << 2) & 63), alphabet);
            case 1:
                unpacked[i] = util::match_letter<PROTO_OUT>((packed(in_byte)) & 63, alphabet);
        }
    }


    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    void unpack_common(const Sequence<INTERNAL_IN> &packed,
                       ProtoSequence<INTERNAL_OUT, PROTO_OUT> &unpacked,
                       const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                unpack_common_2(packed, unpacked, alphabet);
                return;
            case 3:
                unpack_common_3(packed, unpacked, alphabet);
                return;
            case 4:
                unpack_common_4(packed, unpacked, alphabet);
                return;
            case 5:
                unpack_common_5(packed, unpacked, alphabet);
                return;
            case 6:
                unpack_common_6(packed, unpacked, alphabet);
                return;
            default:
                throw std::invalid_argument(std::string("\"alphabet\" has invalid alphabet size - it is ") +
                                                    std::to_string(alphabet.alphabet_size()) +
                                                    " but it should be between 2 and 6 inclusive");
        }
    }
}

#define FETCH(reg_num) \
    v##reg_num = *it_in;\
    ++it_in;
#define ALIGN_0(reg_num, trim) \
    *it_out = v##reg_num & trim##u; \
    ++it_out;
#define ALIGN_R(reg_num, shift, trim) \
    *it_out = (v##reg_num >> shift##u) & trim##u; \
    ++it_out;
#define ALIGN_RL(reg_num_a, shift_a, trim_a, reg_num_b, shift_b, trim_b) \
    *it_out = ((v##reg_num_a >> shift_a##u) & trim_a##u) | ((v##reg_num_b << shift_b##u) & trim_b##u); \
    ++it_out;
namespace tidysq::alt::internal {
    template<typename ITER_CONST_IN, typename ITER_OUT>
    void unpack_octet_5(ITER_CONST_IN &it_in, ITER_OUT &it_out) {
        LetterValue v1, v2;
        FETCH(1)
        FETCH(2)
        ALIGN_0(1, 31)
        ALIGN_RL(1, 5, 7, 2, 3, 31)
        ALIGN_R(2, 2, 31)
        FETCH(1)
        ALIGN_RL(2, 7, 1, 1, 1, 31)
        FETCH(2)
        ALIGN_RL(1, 4, 15, 2, 4, 31)
        ALIGN_R(2, 1, 31)
        FETCH(1)
        ALIGN_RL(2, 6, 3, 1, 2, 31)
        ALIGN_R(1, 3, 31)
    }
}

#undef FETCH
#undef ALIGN_0
#undef ALIGN_R
#undef ALIGN_RL
