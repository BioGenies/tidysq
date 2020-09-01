#ifndef TIDYSQ_SQRCPP_H
#define TIDYSQ_SQRCPP_H

#include <utility>

#include "general.h"
#include "SequenceRCPP.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"
#include "tidysq/util/transition.h"
#include "tidysq/util/alphabet.h"

namespace tidysq {
    template<>
    class Sq<RCPP> {
        Rcpp::List content_;
        Alphabet alphabet_;
        SqType type_;
    public:
        typedef Sequence<RCPP> SequenceType;

        Sq(lensq length, Alphabet alphabet, const SqType &type) :
                content_(Rcpp::List(length)),
                alphabet_(std::move(alphabet)),
                type_(type) {};

        Sq(Alphabet alphabet, const SqType &type) :
                Sq(0, std::move(alphabet), type) {};

        Sq(lensq length, const SqType &type) :
                Sq(length, util::getStandardAlphabet<0>(type), type) {};

        explicit Sq(const SqType &type) :
                Sq(util::getStandardAlphabet<0>(type), type) {};

        Sq(lensq length, Alphabet alphabet) :
                Sq(length, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        explicit Sq(Alphabet alphabet) :
                Sq(std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        explicit Sq(const Rcpp::List &content) :
                content_(content),
                alphabet_(Alphabet(content.attr("alphabet"))),
                type_(util::getSqType<0>(content.attr("class"))) {};

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
            content_.attr("class") = util::getClassStringVector<0>(type_);
            return content_;
        }
    };
}

#endif //TIDYSQ_SQRCPP_H
