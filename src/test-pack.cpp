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
  }
}