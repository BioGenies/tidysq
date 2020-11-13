#ifndef TIDYSQ_SQ_H
#define TIDYSQ_SQ_H

#include "tidysq/types/general.h"
#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"
#include "tidysq/ops/OperationUnpack.h"
#include "tidysq/types/Proxy.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSq;

    template<InternalType INTERNAL>
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

        template<InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
        ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<INTERNAL>, ProtoSq<INTERNAL_OUT, PROTO_OUT>>(*this, ops::OperationUnpack<INTERNAL, INTERNAL_OUT, PROTO_OUT>());
        }

        friend Rcpp::List export_to_R(const Sq<RCPP> &sq);
    };

    template<>
    inline void Sq<RCPP>::push_back(const ElementType &sequence) {
        content_.push_back(sequence.content());
    }
}

#endif //TIDYSQ_SQ_H
