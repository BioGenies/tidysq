#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/Alphabet.h"
#include "tidysq/ops/unpack.h"
#include "tidysq/Proxy.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    template<typename INTERNAL, typename PROTO>
    class ProtoSq;

    template<typename INTERNAL>
    class Sq {
        typename INTERNAL::SqContentStorageType content_;
        Alphabet alphabet_;
    public:
        typedef Sequence<INTERNAL>                      ElementType;
        typedef typename INTERNAL::SqContentStorageType ContentStorageType;

        Sq(const ContentStorageType &content, const Alphabet &alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        Sq(const LenSq length, const Alphabet &alphabet) :
                Sq(ContentStorageType(length), alphabet) {};

        explicit Sq(const Alphabet &alphabet) :
                Sq(0, alphabet) {};

        Sq(const LenSq length, const SqType &type) :
                Sq(length, Alphabet(type)) {};

        explicit Sq(const SqType &type) :
                Sq(0, Alphabet(type)) {};

        inline SequenceProxy<INTERNAL> operator[](const LenSq index) {
            return SequenceProxy<INTERNAL>(content_[index]);
        }

        inline SequenceConstProxy<INTERNAL> operator[](const LenSq index) const {
            return SequenceConstProxy<INTERNAL>(content_[index]);
        }

        [[nodiscard]] LenSq size() const {
            return content_.size();
        }

        [[nodiscard]] inline const Alphabet& alphabet() const {
            return alphabet_;
        }

        [[nodiscard]] inline const SqType& type() const {
            return alphabet_.type();
        }

        inline void push_back(const ElementType &sequence) {
            content_.push_back(sequence);
        }

        inline bool operator==(const Sq<INTERNAL> &other) const {
            if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
            for (u_LenSq i = 0; i < content_.size(); i++) {
                if ((*this)[i] != other[i]) return false;
            }
            return true;
        }

        inline bool operator!=(const Sq<INTERNAL> &other) const {
            return !operator==(other);
        }

        template<typename INTERNAL_OUT, typename PROTO_OUT>
        ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return tidysq::unpack<INTERNAL, INTERNAL_OUT, PROTO_OUT>(*this);
        }

        template<typename INTERNAL_OUT, typename PROTO_OUT>
        ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack(const LenSq from, const LenSq to) const {
            return tidysq::unpack<INTERNAL, INTERNAL_OUT, PROTO_OUT>(*this, from, to);
        }

        friend Rcpp::List export_to_R(const Sq<RCPP_IT> &sq);
    };

    template<>
    inline void Sq<RCPP_IT>::push_back(const ElementType &sequence) {
        Rcpp::RawVector content = sequence.content();
        content.attr("original_length") = sequence.original_length();
        content_.push_back(content);
    }
}

