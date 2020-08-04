#ifndef TIDYSQ_TYPES_RCPPSQ_H
#define TIDYSQ_TYPES_RCPPSQ_H

#include "types_RcppSequence.h"
#include "types_general.h"

namespace tidysq {
    class RcppSq {
        Rcpp::List content_;
        RcppAlphabet alphabet_;
    public:
        typedef RcppSequence SequenceType;
        typedef RcppAlphabet AlphabetType;

        explicit RcppSq(lensq length, AlphabetType alphabet) :
                content_(Rcpp::List(length)),
                alphabet_(alphabet) {};

        RcppSq(Rcpp::List content, AlphabetType alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        Rcpp::List::Proxy operator[] (lensq index) {
            return content_[index];
        }

        Rcpp::List::const_Proxy operator[] (lensq index) const {
            return content_[index];
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }

    };
}

#endif //TIDYSQ_TYPES_RCPPSQ_H
