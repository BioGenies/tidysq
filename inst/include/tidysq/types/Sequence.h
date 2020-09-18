#ifndef TIDYSQ_SEQUENCE_H
#define TIDYSQ_SEQUENCE_H

#include <vector>
#include <Rcpp.h>

#include "tidysq/types/general.h"
#include "tidysq/types/SequenceProto.h"

namespace tidysq {
    template<>
    class Sequence<ANY_INTERNAL> {
        lensq originalLength_;
    protected:
        Sequence() :
                Sequence(0) {};

        explicit Sequence(lensq originalLength) :
                originalLength_(originalLength) {};

    public:
        typedef ElemRaws ElementType;

        Sequence(const Sequence &other) noexcept = default;

        Sequence(Sequence &&other) noexcept = default;

        virtual const ElementType &operator[](lensq index) const = 0;

        virtual ElementType &operator[](lensq index) = 0;

        [[nodiscard]] lensq originalLength() const {
            return originalLength_;
        }

        [[nodiscard]] virtual lensq size() const = 0;
    };

    template<>
    class Sequence<STD> : public Sequence<ANY_INTERNAL> {
        std::vector<unsigned char> content_;
    public:
        typedef unsigned char ElementType;

        Sequence() :
                Sequence(0, 0) {};

        Sequence(lensq packedLength, lensq originalLength) :
                Sequence<ANY_INTERNAL>(originalLength),
                content_(packedLength) {};

        Sequence(const Sequence &other) noexcept = default;

        Sequence(Sequence &&other) noexcept = default;

        inline const ElementType &operator[](lensq index) const override {
            return content_[index];
        }

        inline ElementType &operator[](lensq index) override {
            return content_[index];
        }

        [[nodiscard]] inline lensq size() const override {
            return content_.size();
        }
    };

    template<>
    class Sequence<RCPP> : public Sequence<ANY_INTERNAL> {
        Rcpp::RawVector content_;
    public:
        typedef unsigned char ElementType;

        Sequence() :
                Sequence(0, 0) {};

        Sequence(lensq packedLength, lensq originalLength) :
                Sequence<ANY_INTERNAL>(originalLength),
                content_(packedLength) {};

        explicit Sequence(const Rcpp::RawVector& content) :
                Sequence<ANY_INTERNAL>(content.attr("original_length")),
                content_(content) {
            if (!content.hasAttribute("original_length"))
                throw std::invalid_argument(R"("content" argument in Sequence should have "original_length" attribute!)");
        };

        Sequence(const Sequence &other) noexcept = default;

        Sequence(Sequence &&other) noexcept = default;

        inline Rcpp::RawVector::const_Proxy operator[] (lensq index) const override {
            return content_(index);
        }

        inline Rcpp::RawVector::Proxy operator[] (lensq index) override {
            return content_(index);
        }

        [[nodiscard]] inline lensq size() const override {
            return content_.size();
        }
    };
}

#endif //TIDYSQ_SEQUENCE_H
