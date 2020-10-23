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

template<>
void test_pack_RCPP<STRING>(const std::vector<typename TypeMapper<RCPP, STRING>::ProtoSequenceContentType> &proto, const Alphabet &alphabet) {
    ProtoSq<RCPP, STRING> sq_proto = ProtoSq<RCPP, STRING>(util::convertStringVector(proto), alphabet);
    ProtoSq<RCPP, STRING> repacked = sq_proto.template pack<RCPP>().template unpack<RCPP, STRING>();

    expect_true(repacked == sq_proto);
}

context("test_packing") {
  test_that("pack RCPP RAWS") {
        test_pack_RCPP<RAWS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP<RAWS>({
                                    {0, 1, 2, 3},
                                    {0, 0},
                                    {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                    {7, 0, 7},
                                    {}
                            }, Alphabet(DNA_BSC));

        test_pack_RCPP<INTS>({{0}},
                             Alphabet(DNA_BSC));

        test_pack_RCPP<INTS>({
                                     {0, 1, 2, 3},
                                     {0, 0},
                                     {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                                     {7, 0, 7},
                                     {}
                             }, Alphabet(DNA_BSC));

        test_pack_RCPP<STRINGS>({{"A"}},
                                Alphabet(DNA_BSC));

        test_pack_RCPP<STRINGS>({
                                     {"A", "C", "G", "T"},
                                     Rcpp::StringVector::create("A", "A"),
                                     {"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"},
                                     Rcpp::StringVector::create("!", "A", "!"),
                                     {}
                             }, Alphabet(DNA_BSC));

        test_pack_RCPP<STRING>({
            "ACTG",
            "AA",
            "GGGGGGGGGGGGGG",
            "!A!",
            ""
        }, Alphabet(DNA_BSC));

        test_pack_RCPP<STRING>({
                                       "XDXdxDxd",
                                       "XDXd",
                                       "XdXdXdXd",
                                       ""
                               }, Alphabet(std::vector<Letter>{"XD", "xd", "xD"}, "Xd"));
  }
}