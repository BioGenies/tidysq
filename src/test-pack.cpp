#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/exports.h"

using namespace tidysq;

context("test_packing") {
  test_that("pack STD RAWS 1") {
    Alphabet alphabet({"A", "C", "G", "T", "-"}, "!");

    ProtoSq<STD, RAWS> sq_proto({
        {'A', 'C', 'T', 'G'},
        {'A', 'A'},
        {'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'}
    }, alphabet);

    Sq<STD> sq = sq_proto.pack<STD>();
  }

  test_that("pack RCPP RAWS 1") {
    Rcpp::StringVector alph_content{"A", "C", "G", "T", "-"};
    alph_content.attr("na_letter") = "!";

    Alphabet alphabet(alph_content);

    ProtoSq<RCPP, RAWS> sq_proto = importProtoFromR<RAWS>(
        Rcpp::List::create(
                Rcpp::RawVector{'A', 'C', 'T', 'G'},
                Rcpp::RawVector{'A', 'A'},
                Rcpp::RawVector{'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'}),
                alph_content);

    Sq<RCPP> sq = sq_proto.pack<RCPP>();
  }
}