#ifndef TIDYSQ_SQPROTORCPP_H
#define TIDYSQ_SQPROTORCPP_H

#include "SequenceProtoRCPP.h"
#include "AlphabetRCPP.h"
#include "../sqapply.h"
#include "../impl_Operation.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<RCPP, PROTO> {
        Rcpp::List content_;
        Alphabet<RCPP> alphabet_;
    public:
        typedef SequenceProto<RCPP, PROTO> SequenceType;
        typedef Alphabet<RCPP> AlphabetType;

        SqProto(Rcpp::List content, AlphabetType alphabet) :
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
            return sqapply<SqProto<RCPP, PROTO>, TYPE_OUT, AlphabetType>(*this, ops::OperationPack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQPROTORCPP_H
