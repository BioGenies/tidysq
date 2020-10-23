#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List& x,
                            const Rcpp::StringVector& motifs) {
    const Sq<RCPP> sq = importFromR(x, "!");
    const Alphabet& alph = sq.alphabet();
    // TODO: die again
    // Steps to take:
    // 1. handle ^ and $ signs
    // 2. convert motifs to a list of lists of letters
    //      (if you go for a class, place it within tidysq-internal namespace)
    //      (each element of inner lists is a set of acceptable letters at this position - stored as bits)
    //      (store one motif as linked list of arrays)
    // 3. iterate over sequences and motifs, trying to align them
    //      (probably for each vector in a list within the list use any_of(letter_vec...))
    return true;
}
