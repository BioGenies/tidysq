#ifndef TIDYSQ_SQPROTORCPP_H
#define TIDYSQ_SQPROTORCPP_H

#include <utility>
#include <tidysq/util/alphabet.h>

#include "SequenceProtoRCPP.h"
#include "Alphabet.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<RCPP, PROTO> {
        Rcpp::List content_;
        Alphabet alphabet_;
        SqType type_;
    public:
        typedef SequenceProto<RCPP, PROTO> SequenceType;

        SqProto(const Rcpp::List &content, Alphabet alphabet, const SqType &type) :
                content_(content),
                alphabet_(std::move(alphabet)),
                type_(type) {};

        SqProto(const lensq length, Alphabet alphabet, const SqType &type) :
                SqProto(Rcpp::List(length), std::move(alphabet), type) {};

        SqProto(const Rcpp::List &content, Alphabet alphabet) :
                SqProto(content, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        SqProto(const lensq length, Alphabet alphabet) :
                SqProto(length, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        SqProto(const Rcpp::List &content, const SqType &type) :
                SqProto(content, util::getStandardAlphabet<0>(type), type) {};

        SqProto(const lensq length, const SqType &type) :
                SqProto(length, util::getStandardAlphabet<0>(type), type) {};

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
        SqType type_;
    public:
        typedef SequenceProto<RCPP, STRING> SequenceType;

        SqProto(const Rcpp::StringVector &content, Alphabet alphabet, const SqType &type) :
                content_(content),
                alphabet_(std::move(alphabet)),
                type_(type) {};

        SqProto(const lensq length, const Alphabet& alphabet, const SqType &type) :
                SqProto(Rcpp::StringVector(length), alphabet, type) {};

        SqProto(const Rcpp::StringVector &content, const Alphabet& alphabet) :
                SqProto(content, alphabet, util::guessSqType<0>(alphabet)) {};

        SqProto(const lensq length, const Alphabet& alphabet) :
                SqProto(length, alphabet, util::guessSqType<0>(alphabet)) {};

        SqProto(const Rcpp::StringVector &content, const SqType &type) :
                SqProto(content, util::getStandardAlphabet<0>(type), type) {};

        SqProto(const lensq length, const SqType &type) :
                SqProto(length, util::getStandardAlphabet<0>(type), type) {};

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
