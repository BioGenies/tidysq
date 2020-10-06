#ifndef TIDYSQ_PROTOSEQUENCE_H
#define TIDYSQ_PROTOSEQUENCE_H

#include <utility>

#include "tidysq/types/general.h"
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

        explicit operator ContentType() {
            return content_;
        }
    };
}

#endif //TIDYSQ_PROTOSEQUENCE_H
