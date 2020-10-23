#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/internal/util.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x, const Rcpp::IntegerVector& indices) {
    const Sq<RCPP> sq = importFromR(x, "!");
    Sq<RCPP> ret(sq.length(), sq.alphabet());
    const AlphSize alph_size = sq.alphabet().alphabet_size();
    bool warning_called = false;
    // TODO: replace with NULL once it works
    Rcpp::StringVector NA_warning = "";

    for (LenSq i = 0; i < sq.length(); ++i) {
        const Sequence<RCPP> sequence = sq[i];
        Sequence<RCPP> out_sequence(
                internal::calculatePackedLength(indices.length(), sq.alphabet()),
                indices.length()
        );

        auto index_iter = indices.begin();
        auto sequence_iter = sequence.begin(sq.alphabet());
        auto out_sequence_iter = out_sequence.begin(sq.alphabet());

        while (index_iter != indices.end() || out_sequence_iter != out_sequence.end(sq.alphabet())) {
            ElementPacked element = 0xffu >> (8u - alph_size);
            if (*index_iter <= sequence.originalLength()) {
                element = sequence_iter.access(*index_iter - 1);
            } else if (!warning_called) {
                NA_warning = "some sequences are subsetted with index bigger than length - NA introduced";
                warning_called = true;
            }
            out_sequence_iter.assign(element);

            ++index_iter;
            ++out_sequence_iter;
        }

        ret[i] = out_sequence;
    }
    return Rcpp::List::create(Rcpp::Named("warning") = NA_warning,
                              Rcpp::Named("sq") = ret.exportToR());
}
