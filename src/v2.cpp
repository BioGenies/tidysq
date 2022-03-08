#include <Rcpp.h>

#include "tidysq/internal/pack.h"

using namespace tidysq;

//[[Rcpp::export]]
std::vector<unsigned char> CPP_test_v2() {
    auto input = std::vector<unsigned char>{0, 1, 2, 3, 4, 5, 6, 7};
    auto input_iter = input.cbegin();
    auto output = std::vector<unsigned char> (3, 0);
    auto output_iter = output.begin();

    tidysq::v2::internal::pack_3<8>(input_iter, output_iter);

    return output;
}