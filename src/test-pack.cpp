#include <testthat.h>

#include "tidysq.h"

using namespace tidysq;

template<typename INTERNAL_UNPACKED, typename PROTO, typename INTERNAL_PACKED>
void test_packing_and_unpacking(const ProtoSq<INTERNAL_UNPACKED, PROTO> &proto_sq) {
    auto packed = proto_sq.template pack<INTERNAL_PACKED>();
    auto reunpacked = packed.template unpack<INTERNAL_UNPACKED, PROTO>();
    expect_true(proto_sq == reunpacked);
}

template<typename INTERNAL, typename PROTO>
ProtoSequence<INTERNAL, PROTO>
create_proto_sequence_from_raws(const std::vector<ElementRaws> &raws, const Alphabet &alphabet);

template<>
ProtoSequence<STD_IT, RAWS_PT>
create_proto_sequence_from_raws<STD_IT, RAWS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    return ProtoSequence<STD_IT, RAWS_PT>(raws);
}

template<>
ProtoSequence<RCPP_IT, RAWS_PT>
create_proto_sequence_from_raws<RCPP_IT, RAWS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    return ProtoSequence<RCPP_IT, RAWS_PT>(Rcpp::wrap(raws));
}

template<>
ProtoSequence<STD_IT, INTS_PT>
create_proto_sequence_from_raws<STD_IT, INTS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    std::vector<ElementInts> ints(raws.size());
    std::copy(raws.cbegin(), raws.cend(), ints.begin());
    return ProtoSequence<STD_IT, INTS_PT>(ints);
}

template<>
ProtoSequence<RCPP_IT, INTS_PT>
create_proto_sequence_from_raws<RCPP_IT, INTS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    Rcpp::IntegerVector ints(raws.size());
    std::copy(raws.cbegin(), raws.cend(), ints.begin());
    return ProtoSequence<RCPP_IT, INTS_PT>(ints);
}

template<>
ProtoSequence<STD_IT, STRINGS_PT>
create_proto_sequence_from_raws<STD_IT, STRINGS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    std::vector<ElementStrings> strings(raws.size());
    for (LenSq i = 0; i < raws.size(); i++) {
        strings[i] = util::match_letter<STRINGS_PT>(raws[i], alphabet);
    }
    return ProtoSequence<STD_IT, STRINGS_PT>(strings);
}

template<>
ProtoSequence<RCPP_IT, STRINGS_PT>
create_proto_sequence_from_raws<RCPP_IT, STRINGS_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    Rcpp::StringVector strings(raws.size());
    for (LenSq i = 0; i < raws.size(); i++) {
        strings[i] = util::match_letter<STRINGS_PT>(raws[i], alphabet);
    }
    return ProtoSequence<RCPP_IT, STRINGS_PT>(strings);
}

template<>
ProtoSequence<STD_IT, STRING_PT>
create_proto_sequence_from_raws<STD_IT, STRING_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    std::string seq;
    for (LenSq i = 0; i < raws.size(); i++) {
        seq += util::match_letter_multichar(raws[i], alphabet);
    }
    return ProtoSequence<STD_IT, STRING_PT>(seq);
}

template<>
ProtoSequence<RCPP_IT, STRING_PT>
create_proto_sequence_from_raws<RCPP_IT, STRING_PT>(const std::vector<ElementRaws> &raws, const Alphabet &alphabet) {
    std::string seq;
    for (LenSq i = 0; i < raws.size(); i++) {
        seq += util::match_letter_multichar(raws[i], alphabet);
    }
    return ProtoSequence<RCPP_IT, STRING_PT>(seq);
}

template<typename INTERNAL, typename PROTO>
ProtoSq<INTERNAL, PROTO> create_proto_sq_from_raws(const std::vector<std::vector<ElementRaws>> &raws, const Alphabet &alphabet) {
    ProtoSq<INTERNAL, PROTO> ret(raws.size(), alphabet);
    for (LenSq i = 0; i < raws.size(); i++) {
        ret[i] = create_proto_sequence_from_raws<INTERNAL, PROTO>(raws[i], alphabet);
    }
    return ret;
}

context("test-pack") {
    const std::vector<std::vector<ElementRaws>> alph_2_sequences{
            {0, 1, 2},
            {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {},
            {0, 3, 0, 3}
    };

    const Alphabet alph_2_simple({"A", "B", "C"});
    const Alphabet alph_2_multichar({"AA", "B", "CCC"});

    const std::vector<std::vector<ElementRaws>> alph_3_sequences{
            {0, 1, 2, 3, 4},
            {0, 0},
            {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
            {},
            {0, 7, 0, 7}
    };

    const Alphabet alph_3_simple = Alphabet(DNA_BSC);
    const Alphabet alph_3_multichar({"AA", "B", "CCC", "D", "EE", "FFF", "G"}, "???");

    std::vector<std::vector<ElementRaws>> alph_4_sequences{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
            {0, 0},
            {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
            {},
            {0, 15, 0, 15}
    };

    const Alphabet alph_4_simple({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"});
    const Alphabet alph_4_multichar({"Aa" "Bb" "Cc" "Dd" "Ee" "Ff" "Gg" "Hh" "Ii" "Jj" "Kk" "Ll" "Mm" "Nn" "Oo"});

    const std::vector<std::vector<ElementRaws>> alph_5_sequences{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
             29, 30},
            {0, 0},
            {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
            {},
            {0, 31, 0, 31}
    };

    const Alphabet alph_5_simple({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q",
                                  "r", "s", "t", "u", "v", "w", "x", "y", "z", "A", "B", "C", "D", "E"}, "!", false);
    const Alphabet alph_5_multichar({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q",
                                     "r", "s", "t", "u", "v", "w", "x", "y", "z", "aa", "bb", "cc", "dd", "ee"});

    const std::vector<std::vector<ElementRaws>> alph_6_sequences{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
             29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
             56, 57, 58, 59, 60, 61},
            {0, 0},
            {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
            {},
            {0, 31, 0, 31}
    };

    const Alphabet alph_6_simple({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q",
                                  "r", "s", "t", "u", "v", "w", "x", "y", "z", "A", "B", "C", "D", "E", "F", "G", "H",
                                  "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y",
                                  "Z", "!", "@", "#", "$", "%", "^", "&", "*", "?", "-", "+"}, "!", false);

    const Alphabet alph_6_multichar({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q",
                                  "r", "s", "t", "u", "v", "w", "x", "y", "z", "A", "B", "C", "D", "E", "F", "G", "H",
                                  "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y",
                                  "Z", "!", "@", "#", "$", "%", "^", "&", "*", "?", "-", "+"}, "NA");




#define TEST_PACK_AND_UNPACK_BASIC(INTERNAL_UNPACKED, PROTO, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)                       \
test_that(std::string("pack and unpack ") + #INTERNAL_UNPACKED "_IT, " + #PROTO +"_PT -> " + #INTERNAL_PACKED +        \
          "_IT, alph_size = " + #ALPH_SIZE + ", simple = " + #SIMPLE) {                                                \
    test_packing_and_unpacking<INTERNAL_UNPACKED##_IT, PROTO##_PT, INTERNAL_PACKED##_IT>(                              \
        create_proto_sq_from_raws<INTERNAL_UNPACKED##_IT, PROTO##_PT>(alph_##ALPH_SIZE##_sequences,                    \
        SIMPLE ? alph_##ALPH_SIZE##_simple : alph_##ALPH_SIZE##_multichar)                                             \
    );                                                                                                                 \
}

#define TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)                      \
    TEST_PACK_AND_UNPACK_BASIC(INTERNAL_UNPACKED, RAWS, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)                            \
    TEST_PACK_AND_UNPACK_BASIC(INTERNAL_UNPACKED, INTS, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)                            \
    TEST_PACK_AND_UNPACK_BASIC(INTERNAL_UNPACKED, STRINGS, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)                         \
    TEST_PACK_AND_UNPACK_BASIC(INTERNAL_UNPACKED, STRING, INTERNAL_PACKED, ALPH_SIZE, SIMPLE)

#define TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(INTERNAL_UNPACKED, INTERNAL_PACKED, SIMPLE)              \
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, 2, SIMPLE)                                  \
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, 3, SIMPLE)                                  \
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, 4, SIMPLE)                                  \
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, 5, SIMPLE)                                  \
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO(INTERNAL_UNPACKED, INTERNAL_PACKED, 6, SIMPLE)                                  \

    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(RCPP, RCPP, true)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(RCPP, RCPP, false)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(RCPP, STD, true)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(RCPP, STD, false)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(STD, RCPP, true)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(STD, RCPP, false)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(STD, STD, true)
    TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES(STD, STD, false)


#undef TEST_PACK_AND_UNPACK_BASIC
#undef TEST_PACK_AND_UNPACK_FOR_ALL_PROTO
#undef TEST_PACK_AND_UNPACK_FOR_ALL_PROTO_FOR_ALL_ALPH_SIZES

}
