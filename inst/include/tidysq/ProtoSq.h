#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/Alphabet.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/pack.h"
#include "tidysq/Proxy.h"

namespace tidysq {
    template<typename INTERNAL>
    class Sq;

    template<typename INTERNAL, typename PROTO>
    inline typename ProtoSq<INTERNAL, PROTO>::ContentStorageType export_to_R(const ProtoSq<INTERNAL, PROTO> &proto_sq);

    template<typename INTERNAL, typename PROTO>
    class ProtoSq {
        typename util::TypeBinder<INTERNAL, PROTO>::ProtoSqContentStorageType content_;
        Alphabet alphabet_;
    public:
        typedef ProtoSequence<INTERNAL, PROTO>                                      ElementType;
        typedef typename util::TypeBinder<INTERNAL, PROTO>::ProtoSqContentStorageType     ContentStorageType;

        ProtoSq(const ContentStorageType &content, const Alphabet &alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        ProtoSq(const LenSq length, const Alphabet &alphabet) :
                ProtoSq(ContentStorageType(length), alphabet) {};

        ProtoSq(const ContentStorageType &content, const SqType &type) :
                ProtoSq(content, Alphabet(type)) {};

        ProtoSq(const LenSq length, const SqType &type) :
                ProtoSq(length, Alphabet(type)) {};

        inline ProtoSequenceProxy<INTERNAL, PROTO> operator[](const LenSq index) {
            return ProtoSequenceProxy<INTERNAL, PROTO>(content_[index]);
        }

        inline ProtoSequenceConstProxy<INTERNAL, PROTO> operator[](const LenSq index) const {
            return ProtoSequenceConstProxy<INTERNAL, PROTO>(content_[index]);
        }

        [[nodiscard]] inline LenSq size() const {
            return content_.size();
        }

        [[nodiscard]] inline const Alphabet &alphabet() const {
            return alphabet_;
        }

        inline Alphabet& alphabet() {
            return alphabet_;
        }

        [[nodiscard]] inline const SqType &type() const {
            return alphabet_.type();
        }

        template<typename INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return tidysq::pack<INTERNAL, PROTO, INTERNAL_OUT>(*this);
        }

        inline bool operator==(const ProtoSq<INTERNAL, PROTO> &other) const {
           if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
           for (u_LenSq i = 0; i < content_.size(); i++) {
               if ((*this)[i] != other[i]) return false;
           }
           return true;
        }

        inline bool operator!=(const ProtoSq<INTERNAL, PROTO> &other) const {
            return !operator==(other);
        }


        friend typename ProtoSq<INTERNAL, PROTO>::ContentStorageType export_to_R<INTERNAL, PROTO>(const ProtoSq<INTERNAL, PROTO> &proto_sq);
    };


    template<>
    inline bool ProtoSq<RCPP_IT, STRING_PT>::operator==(const ProtoSq<RCPP_IT, STRING_PT> &other) const {
        if ((alphabet_ != other.alphabet_) || (content_.size() != other.content_.size())) return false;
        for (u_LenSq i = 0; i < content_.size(); i++) {
            if ((*this)[i] != other[i]) return false;
        }
        return true;
    }


}
