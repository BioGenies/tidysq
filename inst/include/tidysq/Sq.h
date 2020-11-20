#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/Alphabet.h"
#include "tidysq/TypeMapper.h"
#include "tidysq/ops/OperationUnpack.h"
#include "tidysq/Proxy.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    template<typename INTERNAL, typename PROTO>
    class ProtoSq;

    template<typename INTERNAL>
    class Sq {
        typename InternalTypeMapper<INTERNAL>::SqContentType content_;
        Alphabet alphabet_;
    public:
        typedef typename InternalTypeMapper<INTERNAL>::SqContentType ContentType;
        typedef typename InternalTypeMapper<INTERNAL>::SqElementType ElementType;

        Sq(const ContentType &content, const Alphabet &alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        Sq(const LenSq length, const Alphabet &alphabet) :
                Sq(ContentType(length), alphabet) {};

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

        [[nodiscard]] LenSq length() const {
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

        inline bool operator==(const Sq<INTERNAL> &other) {
            if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
            for (LenSq i = 0; i < content_.size(); i++) {
                if ((*this)[i] != other[i]) return false;
            }
            return true;
        }

        inline bool operator!=(const Sq<INTERNAL> &other) {
            return !operator==(other);
        }

        template<typename INTERNAL_OUT, typename PROTO_OUT>
        ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<INTERNAL>, ProtoSq<INTERNAL_OUT, PROTO_OUT>>(*this, ops::OperationUnpack<INTERNAL, INTERNAL_OUT, PROTO_OUT>());
        }

        friend Rcpp::List export_to_R(const Sq<RCPP_IT> &sq);
    };

    template<>
    inline void Sq<RCPP_IT>::push_back(const ElementType &sequence) {
        content_.push_back(sequence.content());
    }
}

