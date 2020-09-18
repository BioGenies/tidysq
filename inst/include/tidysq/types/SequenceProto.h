#ifndef TIDYSQ_SEQUENCEPROTO_H
#define TIDYSQ_SEQUENCEPROTO_H

#include <utility>
#include <vector>
#include "general.h"

namespace tidysq {
    template<>
    class SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    protected:
        SequenceProto() noexcept = default;

    public:
        SequenceProto(const SequenceProto &other) noexcept = default;

        SequenceProto(SequenceProto &&other) noexcept = default;

        [[nodiscard]] virtual lensq size() const = 0;

        [[nodiscard]] virtual letvalue getLetterValue(lensq i, const Alphabet &alphabet) const = 0;
    };

    template<>
    class SequenceProto<STD, RAWS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
        std::vector<ElemRaws> content_;
    public:
        typedef ElemRaws ElementType;

        explicit SequenceProto(std::vector<ElementType> content) :
                SequenceProto<ANY_INTERNAL, ANY_PROTO>(),
                content_(std::move(content)) {};

        SequenceProto(std::initializer_list<ElementType> list) :
                SequenceProto<ANY_INTERNAL, ANY_PROTO>(),
                content_(std::vector<ElementType>(list)) {};

        SequenceProto(const SequenceProto &other) noexcept = default;

        SequenceProto(SequenceProto &&other) noexcept = default;

        inline const ElementType &operator[](lensq index) const {
            return content_[index];
        }

        inline ElementType &operator[](lensq index) {
            return content_[index];
        }

        [[nodiscard]] inline lensq size() const override {
            return content_.size();
        }

        [[nodiscard]] inline letvalue getLetterValue(const lensq index, const Alphabet &alphabet) const override {
            return content_[index];
        }
    };

    template<>
    class SequenceProto<STD, INTS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<STD, STRINGS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<STD, STRING> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<RCPP, RAWS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<RCPP, INTS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<RCPP, STRINGS> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };

    template<>
    class SequenceProto<RCPP, STRING> : public SequenceProto<ANY_INTERNAL, ANY_PROTO> {
    };
}

#endif //TIDYSQ_SEQUENCEPROTO_H
