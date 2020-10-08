#ifndef TIDYSQ_PROTOSEQUENCE_H
#define TIDYSQ_PROTOSEQUENCE_H

#include "tidysq/types/general.h"

namespace tidysq {
    template<InternalType, ProtoType>
    class ProtoSequence;
}

#include <RcppCommon.h>

namespace Rcpp {
    template<>
    SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::RAWS> &);
}

#include <utility>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"

namespace tidysq {
    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSequence {
        typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType content_;
    public:
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType ContentType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceElementType ElementType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceAccessType AccessType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceConstAccessType ConstAccessType;

        explicit ProtoSequence(const ContentType &content) :
                content_(content) {};

        explicit ProtoSequence(const LenSq length) :
                ProtoSequence(ContentType(length)) {};

        ProtoSequence() :
                ProtoSequence(0) {};

        ProtoSequence(const std::initializer_list<ElementType> &list) :
                content_(ContentType(list)) {};

        ProtoSequence(const ProtoSequence &other) noexcept = default;

        ProtoSequence(ProtoSequence &&other) noexcept = default;

        ProtoSequence& operator=(const ProtoSequence &other) noexcept = default;

        ProtoSequence& operator=(ProtoSequence &&other) noexcept = default;

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        };

        inline ConstAccessType operator[](const LenSq index) const {
            return content_[index];
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline LetValue getLetterValue(const LetValue index, const Alphabet &alphabet) const {
            return content_[index];
        }

        [[nodiscard]] inline ContentType getContent() const {
            return content_;
        }

        inline bool operator==(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return content_ == other.content_;
        }

        inline bool operator!=(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return !operator==(other);
        }
    };
}

namespace Rcpp {
    template<>
    inline SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::RAWS> &obj) {
        return obj.getContent();
    }
}

#endif //TIDYSQ_PROTOSEQUENCE_H
