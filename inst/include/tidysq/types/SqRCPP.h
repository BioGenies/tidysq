#ifndef TIDYSQ_SQRCPP_H
#define TIDYSQ_SQRCPP_H

#include <utility>

#include "general.h"
#include "SequenceRCPP.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"

namespace tidysq {
    template<>
    class Sq<RCPP> {
        Rcpp::List content_;
        Alphabet alphabet_;
    public:
        typedef Sequence<RCPP> SequenceType;

        explicit Sq(Alphabet alphabet) :
                Sq(0, std::move(alphabet)) {};

        Sq(lensq length, Alphabet alphabet) :
                content_(Rcpp::List(length)),
                alphabet_(std::move(alphabet)) {};

        explicit Sq(const Rcpp::List &content) :
                content_(content),
                alphabet_(Alphabet(content.attr("alphabet"))) {};

        Rcpp::List::Proxy operator[] (lensq index) {
            return content_[index];
        }

        Rcpp::List::const_Proxy operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] const Alphabet& alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT,
                ProtoType PROTO_OUT>
        SqProto<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<RCPP>, SqProto<INTERNAL_OUT, PROTO_OUT>>(*this, ops::OperationUnpack<RCPP, INTERNAL_OUT, PROTO_OUT>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = static_cast<Rcpp::StringVector>(alphabet_);
            return content_;
        }
    };
}

#endif //TIDYSQ_SQRCPP_H
