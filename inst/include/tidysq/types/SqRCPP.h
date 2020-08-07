#ifndef TIDYSQ_SQRCPP_H
#define TIDYSQ_SQRCPP_H

#include <utility>

#include "general.h"
#include "SequenceRCPP.h"
#include "AlphabetRCPP.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"

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
                alphabet_(std::move(alphabet)) {};

        explicit Sq(const Rcpp::List &content) :
                content_(content),
                alphabet_(Alphabet<RCPP>(content.attr("alphabet"))) {};

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

        template<InternalType INTERNAL_OUT,
                ProtoType PROTO_OUT>
        SqProto<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<RCPP>, SqProto<INTERNAL_OUT, PROTO_OUT>, AlphabetType>(*this, ops::OperationUnpack<RCPP, INTERNAL_OUT, PROTO_OUT>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQRCPP_H
