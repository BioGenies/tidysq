#ifndef TIDYSQ_PROTOSQ_H
#define TIDYSQ_PROTOSQ_H

#include "tidysq/types/general.h"
#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"
#include "tidysq/util/alphabet.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/OperationPack.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sq;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSq {
        typename TypeMapper<INTERNAL, PROTO>::ProtoSqContentType content_;
        Alphabet alphabet_;
        SqType type_;
    public:
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqContentType ContentType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqElementType ElementType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqAccessType AccessType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSqConstAccessType ConstAccessType;

        ProtoSq(const ContentType &content, const Alphabet &alphabet, const SqType &type) :
                content_(content),
                alphabet_(alphabet),
                type_(type) {};

        ProtoSq(const LenSq length, const Alphabet &alphabet, const SqType &type) :
                ProtoSq(ContentType(length), alphabet, type) {};

        ProtoSq(const ContentType &content, const Alphabet &alphabet) :
                ProtoSq(content, alphabet, util::guessSqType(alphabet)) {};

        ProtoSq(const LenSq length, const Alphabet &alphabet) :
                ProtoSq(length, alphabet, util::guessSqType(alphabet)) {};

        ProtoSq(const ContentType &content, const SqType &type) :
                ProtoSq(length, util::getStandardAlphabet(type), type) {};

        ProtoSq(const LenSq length, const SqType &type) :
                ProtoSq(length, util::getStandardAlphabet(type), type) {};

        inline AccessType operator[] (const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator[] (const LenSq index) const {
            return content_[index];
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline const Alphabet &alphabet() const {
            return alphabet_;
        }

        [[nodiscard]] inline const SqType &type() const {
            return type_;
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() {
            return sqapply<ProtoSq<INTERNAL, PROTO>, Sq<INTERNAL_OUT>>(*this, ops::OperationPack<INTERNAL, PROTO, INTERNAL_OUT>());
        }
    };

    template<ProtoType PROTO>
    ProtoSq<RCPP, PROTO> importProtoFromR(const Rcpp::List &proto,
                                          const Rcpp::StringVector &alphabet) {
        return ProtoSq<RCPP, PROTO>(proto, Alphabet(alphabet));
    }
}

#endif //TIDYSQ_PROTOSQ_H
