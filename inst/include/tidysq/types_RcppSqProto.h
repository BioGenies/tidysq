#ifndef TIDYSQ_TYPES_RCPPSQPROTO_H
#define TIDYSQ_TYPES_RCPPSQPROTO_H

#include "types_RcppSequence.h"
#include "sqapply.h"
#include "impl_Operation.h"

namespace tidysq {
    template<typename RCPP_SEQUENCE_PROTO>
    class RcppSqProto {
        Rcpp::List content_;
        RcppAlphabet alphabet_;
    public:
        typedef RCPP_SEQUENCE_PROTO SequenceType;
        typedef RcppAlphabet AlphabetType;

        RcppSqProto(Rcpp::List content, RcppAlphabet alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        Rcpp::List::Proxy operator[] (lensq index) {
            return content_[index];
        }

        Rcpp::List::const_Proxy operator[] (lensq index) const {
            return content_[index];
        }

        lensq length() const {
            return content_.size();
        }

        AlphabetType alphabet() const {
            return alphabet_;
        }

        template<typename TYPE_OUT>
        TYPE_OUT pack() const {
            return sqapply<RcppSqProto, TYPE_OUT, AlphabetType>(*this, ops::OperationPack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_TYPES_RCPPSQPROTO_H
