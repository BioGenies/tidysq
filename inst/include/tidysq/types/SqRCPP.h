#ifndef TIDYSQ_SQRCPP_H
#define TIDYSQ_SQRCPP_H

#include "general.h"
#include "SequenceRCPP.h"
#include "AlphabetRCPP.h"
#include "../sqapply.h"
#include "../impl_Operation.h"

namespace tidysq {
    template<>
    class Sq<RCPP> {
        Rcpp::List content_;
        Alphabet<RCPP> alphabet_;
    public:
        typedef Sequence<RCPP> SequenceType;
        typedef Alphabet<RCPP> AlphabetType;

        Sq(lensq length, AlphabetType alphabet) :
                content_(Rcpp::List(length)),
                alphabet_(alphabet) {};

        explicit Sq(Rcpp::List content) :
                content_(content),
                alphabet_(content.attr("alphabet")) {};

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
        TYPE_OUT unpack() const {
            return sqapply<Sq<RCPP>, TYPE_OUT, AlphabetType>(*this, ops::OperationUnpack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQRCPP_H
