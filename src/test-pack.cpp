#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/tidysq-includes.h"

using namespace tidysq;

template<ProtoType PROTO>
void test_pack_RCPP(const std::vector<typename TypeMapper<RCPP, PROTO>::ProtoSequenceContentType> &proto, const Alphabet &alphabet) {
    auto list = Rcpp::List::create();

    for (const auto &sequence : proto) {
        list.push_back(sequence);
    }

    ProtoSq<RCPP, PROTO> sq_proto = ProtoSq<RCPP, PROTO>(list, alphabet);
    ProtoSq<RCPP, PROTO> repacked = sq_proto.template pack<RCPP>().template unpack<RCPP, PROTO>();

            expect_true(repacked == sq_proto);
}

template<>
void test_pack_RCPP<STRING>(const std::vector<typename TypeMapper<RCPP, STRING>::ProtoSequenceContentType> &proto, const Alphabet &alphabet) {
    ProtoSq<RCPP, STRING> sq_proto = ProtoSq<RCPP, STRING>(util::convert_string_vector(proto), alphabet);
    ProtoSq<RCPP, STRING> repacked = sq_proto.template pack<RCPP>().template unpack<RCPP, STRING>();

            expect_true(repacked == sq_proto);
}

template<ProtoType PROTO>
void test_pack_STD(const typename TypeMapper<STD, PROTO>::ProtoSqContentType &proto, const Alphabet &alphabet) {
    ProtoSq<STD, PROTO> sq_proto = ProtoSq<STD, PROTO>(proto, alphabet);
    ProtoSq<STD, PROTO> repacked = sq_proto.template pack<STD>().template unpack<STD, PROTO>();

            expect_true(repacked == sq_proto);
}

context("test_packing") {
            test_that("packing RCPP RAWS") {
        test_pack_RCPP<RAWS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP<RAWS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP INTS") {
        test_pack_RCPP<INTS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP<INTS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP STRINGS") {
        test_pack_RCPP<STRINGS>({{"A"}},
                                Alphabet(DNA_BSC));

        test_pack_RCPP<STRINGS>({
                                        {"A", "C", "G", "T"},
                                        Rcpp::StringVector::create("A", "A"),
                                        {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                        Rcpp::StringVector::create("!", "A", "!"),
                                        {}
                                }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP STRING") {
        test_pack_RCPP<STRING>({{"A"}},
                               Alphabet(DNA_BSC));

        test_pack_RCPP<STRING>({
                                       "ACTG",
                                       "AA",
                                       "GGGGGGGGGGGGGG",
                                       "!A!",
                                       ""
                               }, Alphabet(DNA_BSC));
    }
            test_that("packing RCPP STRING MULTICHAR") {
        test_pack_RCPP<STRING>({
                                       "AAmAALJmAmAAmA",
                                       "AJ?mA?J",
                                       "mA",
                                       ""
                               }, Alphabet(std::vector<Letter>{"A", "mA", "L", "J"}, "?"));
    }

            test_that("packing STD RAWS") {
        test_pack_STD<RAWS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_STD<RAWS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing STD INTS") {
        test_pack_STD<INTS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_STD<INTS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));
    }
            test_that("packing STD STRINGS") {
        test_pack_STD<STRINGS>({{"A"}},
                                Alphabet(DNA_BSC));

        test_pack_STD<STRINGS>({
                                        {"A", "C", "G", "T"},
                                        {"A", "A"},
                                        {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                        {"!", "A", "!"},
                                        {}
                                }, Alphabet(DNA_BSC));
    }
            test_that("packing STD STRING") {
        test_pack_STD<STRING>({ProtoSequence<STD, STRING>("A")},
                              Alphabet(DNA_BSC));

        test_pack_STD<STRING>({ProtoSequence<STD, STRING>("ACTG"),
                               ProtoSequence<STD, STRING>("AA"),
                               ProtoSequence<STD, STRING>("GGGGGGGGGGGGGG"),
                               ProtoSequence<STD, STRING>("!A!"),
                               ProtoSequence<STD, STRING>("")
                              }, Alphabet(DNA_BSC));
    }
            test_that("packing STD STRING MULTICHAR") {
        test_pack_STD<STRING>({ProtoSequence<STD, STRING>("AAmAALJmAmAAmA"),
                               ProtoSequence<STD, STRING>("AJ?mA?J"),
                               ProtoSequence<STD, STRING>("mA"),
                               ProtoSequence<STD, STRING>("")
                              }, Alphabet(std::vector<Letter>{"A", "mA", "L", "J"}, "?"));
    }

}