#ifndef TIDYSQ_TYPES_RCPPSEQUENCE_H
#define TIDYSQ_TYPES_RCPPSEQUENCE_H

#include <Rcpp.h>

namespace tidysq {
    typedef Rcpp::RawVector RcppSequence;

    typedef Rcpp::RawVector RcppSequenceProtoRaw;
    typedef Rcpp::IntegerVector RcppSequenceProtoInteger;
    typedef Rcpp::StringVector RcppSequenceProtoString;
    typedef Rcpp::StringVector RcppSequenceProtoStrings;

    typedef Rcpp::StringVector RcppAlphabet;
}

#endif //TIDYSQ_TYPES_RCPPSEQUENCE_H
