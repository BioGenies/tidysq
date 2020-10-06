#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/exports.h"

using namespace tidysq;

context("test_packing") {
//  test_that("pack STD RAWS 1") {
//    ProtoSequence<STD, RAWS> sequence_proto_1{'A', 'C', 'T', 'G'};
//    ProtoSequence<STD, RAWS> sequence_proto_2{'A', 'A'};
//    ProtoSequence<STD, RAWS> sequence_proto_3{'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'};
//
//    Alphabet alphabet(std::vector<std::string>{"A", "C", "G", "T", "-"}, "!");
//
//    ProtoSq<STD, RAWS> sq_proto(
//        std::vector<ProtoSequence<STD, RAWS>>{sequence_proto_1,
//                                              sequence_proto_2,
//                                              sequence_proto_3},
//        alphabet);
//
//    Sq<STD> sq = sq_proto.pack<STD>();
//
//  }
//
//  test_that("pack RCPP RAWS 1") {
//    ProtoSequence<RCPP, RAWS> sequence_proto_1{'A', 'C', 'T', 'G'};
//    ProtoSequence<RCPP, RAWS> sequence_proto_2{'A', 'A'};
//    ProtoSequence<RCPP, RAWS> sequence_proto_3{'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'};
//
//    Rcpp::StringVector alph_content{"A", "C", "G", "T", "-"};
//    alph_content.attr("na_letter") = "!";
//
//    Alphabet alphabet(alph_content);
//
//    ProtoSq<RCPP, RAWS> sq_proto = importProtoFromR<RAWS>(
//        Rcpp::List{sequence_proto_1,
//                   sequence_proto_2,
//                   sequence_proto_3},
//        alph_content);
//
//    Sq<RCPP> sq = sq_proto.pack<RCPP>();
//  }
}