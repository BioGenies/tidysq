#ifndef TIDYSQ_SQPROTORCPP_H
#define TIDYSQ_SQPROTORCPP_H

#include <utility>

#include "SequenceProtoRCPP.h"
#include "Alphabet.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<RCPP, PROTO> {
        Rcpp::List content_;
        Alphabet alphabet_;
    public:
        typedef SequenceProto<RCPP, PROTO> SequenceType;

        SqProto(const Rcpp::List& content, Alphabet alphabet) :
                content_(content),
                alphabet_(std::move(alphabet)) {};

        SqProto(lensq length, Alphabet alphabet) :
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

        [[nodiscard]] const Alphabet& alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<RCPP, PROTO>, Sq<INTERNAL_OUT>>(*this, ops::OperationPack<RCPP, PROTO, INTERNAL_OUT>());
        }

        Rcpp::List exportToR() {
            return content_;
        }
    };

    template<>
    class SqProto<RCPP, STRING> {
        Rcpp::StringVector content_;
        Alphabet alphabet_;
    public:
        typedef SequenceProto<RCPP, STRING> SequenceType;

        SqProto(const Rcpp::StringVector& content, Alphabet alphabet) :
                content_(content),
                alphabet_(std::move(alphabet)) {};

        SqProto(lensq length, Alphabet alphabet) :
                content_(Rcpp::StringVector(length)),
                alphabet_(std::move(alphabet)) {};

        Rcpp::StringVector::Proxy operator[] (lensq index) {
            return content_[index];
        }

        Rcpp::StringVector::const_Proxy operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] const Alphabet &alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<RCPP, STRING>, Sq<INTERNAL_OUT>>(*this, ops::OperationPack<RCPP, STRING, INTERNAL_OUT>());
        }

        Rcpp::StringVector exportToR() {
            return content_;
        }
    };
}

#endif //TIDYSQ_SQPROTORCPP_H
