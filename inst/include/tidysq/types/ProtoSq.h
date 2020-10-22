#ifndef TIDYSQ_PROTOSQ_H
#define TIDYSQ_PROTOSQ_H

#include "tidysq/types/general.h"
#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/OperationPack.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sq;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSq {
        typename TypeMapper<INTERNAL, PROTO>::ProtoSqContentType content_;
        Alphabet alphabet_;
    public:
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqContentType ContentType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqElementType ElementType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType ElementUnderlyingType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqAccessType AccessType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqConstAccessType ConstAccessType;

        ProtoSq(const ContentType &content, const Alphabet &alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        ProtoSq(const LenSq length, const Alphabet &alphabet) :
                ProtoSq(ContentType(length), alphabet) {};

        ProtoSq(const ContentType &content, const SqType &type) :
                ProtoSq(length, Alphabet(type)) {};

        ProtoSq(const LenSq length, const SqType &type) :
                ProtoSq(length, Alphabet(type)) {};

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator[](const LenSq index) const {
            return content_[index];
        }

        inline ElementType get(const LenSq index) const {
            return content_[index];
        }

        inline void set(const LenSq index, const ElementType &value) {
            content_[index] = value;
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline const Alphabet &alphabet() const {
            return alphabet_;
        }

        [[nodiscard]] inline const SqType &type() const {
            return alphabet_.type();
        }

        inline ContentType exportToR() {
            throw std::exception();
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() {
            return sqapply<ProtoSq<INTERNAL, PROTO>, Sq<INTERNAL_OUT>>(*this,
                                                                       ops::OperationPack<INTERNAL, PROTO, INTERNAL_OUT>());
        }

        inline bool operator==(const ProtoSq<INTERNAL, PROTO> &other) {
           if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
           for (LenSq i = 0; i < content_.size(); i++) {
               if (!Rcpp::is_true(Rcpp::all(ElementUnderlyingType(content_[i]) == ElementUnderlyingType(other.content_[i])))) return false;
           }
           return true;
        }

        inline bool operator!=(const ProtoSq<INTERNAL, PROTO> &other) {
            return !operator==(other);
        }
    };

    template<>
    inline ProtoSq<RCPP, RAWS>::ContentType ProtoSq<RCPP, RAWS>::exportToR() {
        return content_;
    }

    template<>
    inline ProtoSq<RCPP, INTS>::ContentType ProtoSq<RCPP, INTS>::exportToR() {
        return content_;
    }

    template<>
    inline ProtoSq<RCPP, STRINGS>::ContentType ProtoSq<RCPP, STRINGS>::exportToR() {
        return content_;
    }

    template<>
    inline ProtoSq<RCPP, STRING>::ContentType ProtoSq<RCPP, STRING>::exportToR() {
        return content_;
    }

    template<>
    inline ProtoSq<RCPP, STRING>::ElementType ProtoSq<RCPP, STRING>::get(const LenSq index) const {
        return ProtoSq<RCPP, STRING>::ElementType(static_cast<std::string>(Rcpp::StringVector(content_[index].get())[0]));
    }

    template<>
    inline void ProtoSq<RCPP, STRING>::set(const LenSq index, const ElementType &value) {
        content_[index] = value.content();
    }

    template<>
    inline bool ProtoSq<RCPP, STRING>::operator==(const ProtoSq<RCPP, STRING> &other) {
        if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
        for (LenSq i = 0; i < content_.size(); i++) {
            if (ElementUnderlyingType(content_[i]) != ElementUnderlyingType(other.content_[i])) return false;
        }
        return true;
    }

    template<ProtoType PROTO>
    ProtoSq<RCPP, PROTO> importProtoFromR(const typename ProtoSq<RCPP, PROTO>::ContentType &proto,
                                          const Rcpp::StringVector &alphabet,
                                          const Rcpp::StringVector &NA_letter) {
        return ProtoSq<RCPP, PROTO>(proto, Alphabet(alphabet, NA_letter));
    }
}

#endif //TIDYSQ_PROTOSQ_H
