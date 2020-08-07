#ifndef TIDYSQ_SQPROTORCPP_H
#define TIDYSQ_SQPROTORCPP_H

#include <utility>

#include "SequenceProtoRCPP.h"
#include "AlphabetRCPP.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<RCPP, PROTO> {
        Rcpp::List content_;
        Alphabet<RCPP> alphabet_;
    public:
        typedef SequenceProto<RCPP, PROTO> SequenceType;
        typedef Alphabet<RCPP> AlphabetType;

        SqProto(const Rcpp::List& content, AlphabetType alphabet) :
                content_(content),
                alphabet_(std::move(alphabet)) {};

        SqProto(lensq length, AlphabetType alphabet) :
                content_(Rcpp::List(length)),
                alphabet_(std::move(alphabet)) {};

        Rcpp::List::Proxy operator[] (lensq index) {
            return content_[index];
        }

        Rcpp::List::const_Proxy operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] AlphabetType alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<RCPP, PROTO>, Sq<INTERNAL_OUT>, AlphabetType>(*this, ops::OperationPack<RCPP, PROTO, INTERNAL_OUT>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQPROTORCPP_H
