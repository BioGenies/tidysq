#include <Rcpp.h>
#include <testthat.h>

#include "tidysq.h"

using namespace tidysq;

template<typename PROTO>
void test_pack_RCPP_IT(const std::vector<typename TypeBinder<RCPP_IT, PROTO>::ProtoSequenceContentStorageType> &proto, const Alphabet &alphabet) {
    auto list = Rcpp::List::create();

    for (const auto &sequence : proto) {
        list.push_back(sequence);
    }

    ProtoSq<RCPP_IT, PROTO> sq_proto = ProtoSq<RCPP_IT, PROTO>(list, alphabet);
    ProtoSq<RCPP_IT, PROTO> repacked = sq_proto.template pack<RCPP_IT>().template unpack<RCPP_IT, PROTO>();

            expect_true(repacked == sq_proto);
}

template<>
void test_pack_RCPP_IT<STRING_PT>(const std::vector<typename TypeBinder<RCPP_IT, STRING_PT>::ProtoSequenceContentStorageType> &proto, const Alphabet &alphabet) {
    ProtoSq<RCPP_IT, STRING_PT> sq_proto = ProtoSq<RCPP_IT, STRING_PT>(util::convert_string_vector(proto), alphabet);
    ProtoSq<RCPP_IT, STRING_PT> repacked = sq_proto.template pack<RCPP_IT>().template unpack<RCPP_IT, STRING_PT>();

            expect_true(repacked == sq_proto);
}

template<typename PROTO>
void test_pack_STD_IT(const typename TypeBinder<STD_IT, PROTO>::ProtoSqContentStorageType &proto, const Alphabet &alphabet) {
    ProtoSq<STD_IT, PROTO> sq_proto = ProtoSq<STD_IT, PROTO>(proto, alphabet);
    ProtoSq<STD_IT, PROTO> repacked = sq_proto.template pack<STD_IT>().template unpack<STD_IT, PROTO>();

            expect_true(repacked == sq_proto);
}

context("test_packing") {
            test_that("packing RCPP_IT RAWS") {
        test_pack_RCPP_IT<RAWS_PT>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP_IT<RAWS_PT>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP_IT INTS") {
        test_pack_RCPP_IT<INTS_PT>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP_IT<INTS_PT>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP_IT STRINGS") {
        test_pack_RCPP_IT<STRINGS_PT>({{"A"}},
                                Alphabet(DNA_BSC));

        test_pack_RCPP_IT<STRINGS_PT>({
                                        {"A", "C", "G", "T"},
                                        Rcpp::StringVector::create("A", "A"),
                                        {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                        Rcpp::StringVector::create("!", "A", "!"),
                                        {}
                                }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP_IT STRING") {
        test_pack_RCPP_IT<STRING_PT>({{"A"}},
                               Alphabet(DNA_BSC));

        test_pack_RCPP_IT<STRING_PT>({
                                       "ACTG",
                                       "AA",
                                       "GGGGGGGGGGGGGG",
                                       "!A!",
                                       ""
                               }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP_IT STRING MULTICHAR") {
        test_pack_RCPP_IT<STRING_PT>({
                                       "AAmAALJmAmAAmA",
                                       "AJ?mA?J",
                                       "mA",
                                       ""
                               }, Alphabet(std::vector<Letter>{"A", "mA", "L", "J"}, "?"));
    }

            test_that("packing STD_IT RAWS") {
        test_pack_STD_IT<RAWS_PT>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_STD_IT<RAWS_PT>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing STD_IT INTS") {
        test_pack_STD_IT<INTS_PT>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_STD_IT<INTS_PT>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing STD_IT STRINGS") {
        test_pack_STD_IT<STRINGS_PT>({{"A"}},
                                Alphabet(DNA_BSC));

        test_pack_STD_IT<STRINGS_PT>({
                                        {"A", "C", "G", "T"},
                                        {"A", "A"},
                                        {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                        {"!", "A", "!"},
                                        {}
                                }, Alphabet(DNA_BSC));
    }
            test_that("packing STD_IT STRING") {
        test_pack_STD_IT<STRING_PT>({ProtoSequence<STD_IT, STRING_PT>("A")},
                              Alphabet(DNA_BSC));

        test_pack_STD_IT<STRING_PT>({ProtoSequence<STD_IT, STRING_PT>("ACTG"),
                               ProtoSequence<STD_IT, STRING_PT>("AA"),
                               ProtoSequence<STD_IT, STRING_PT>("GGGGGGGGGGGGGG"),
                               ProtoSequence<STD_IT, STRING_PT>("!A!"),
                               ProtoSequence<STD_IT, STRING_PT>("")
                              }, Alphabet(DNA_BSC));
    }
            test_that("packing STD_IT STRING MULTICHAR") {
        test_pack_STD_IT<STRING_PT>({ProtoSequence<STD_IT, STRING_PT>("AAmAALJmAmAAmA"),
                               ProtoSequence<STD_IT, STRING_PT>("AJ?mA?J"),
                               ProtoSequence<STD_IT, STRING_PT>("mA"),
                               ProtoSequence<STD_IT, STRING_PT>("")
                              }, Alphabet(std::vector<Letter>{"A", "mA", "L", "J"}, "?"));
    }

}