#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/exports.h"

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

context("test_packing") {
  test_that("pack RCPP RAWS") {
        test_pack_RCPP<RAWS>({{0}},
                            util::getStandardAlphabet(DNA_CLN));

        test_pack_RCPP<RAWS>({
                                    {0, 1, 2, 3},
                                    {0, 0},
                                    {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                    {7, 0, 7},
                                    {}
                            }, util::getStandardAlphabet(DNA_CLN));

        test_pack_RCPP<INTS>({{0}},
                             util::getStandardAlphabet(DNA_CLN));

        test_pack_RCPP<INTS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, util::getStandardAlphabet(DNA_CLN));

        test_pack_RCPP<STRINGS>({{"A"}},
                                util::getStandardAlphabet(DNA_CLN));

        test_pack_RCPP<STRINGS>({
                                     {"A", "C", "G", "T"},
                                     Rcpp::StringVector::create("A", "A"),
                                     {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                     Rcpp::StringVector::create("!", "A", "!"),
                                     {}
                             }, util::getStandardAlphabet(DNA_CLN));

//        test_pack_RCPP<STRING>({
//            "ACTG",
//            "AA",
//            "GGGGGGGGGGGGGG",
//            "!A!",
//            ""
//        }, util::getStandardAlphabet(DNA_CLN));
  }
}