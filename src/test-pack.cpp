#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/exports.h"

using namespace tidysq;

void test_packing_raws(const std::vector<Rcpp::RawVector> &raws, const Alphabet &alphabet) {
    auto list = Rcpp::List::create();

    for (const auto & raw : raws) {
        list.push_back(raw);
    }

    ProtoSq<RCPP, RAWS> sq_proto = ProtoSq<RCPP, RAWS>(list, alphabet);
    ProtoSq<RCPP, RAWS> repacked = sq_proto.pack<RCPP>().unpack<RCPP, RAWS>();

    expect_true(repacked == sq_proto);

}

context("test_packing") {
  test_that("pack RCPP RAWS") {
      test_packing_raws({{1}},
                        util::getStandardAlphabet(DNA_CLN));

      test_packing_raws({
          {1, 2, 3, 4},
          {1, 1},
          {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
          {7, 1, 7},
          {}
      }, util::getStandardAlphabet(DNA_CLN));
  }
}