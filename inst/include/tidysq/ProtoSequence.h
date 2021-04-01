#pragma once

#include <utility>

#include "tidysq/Alphabet.h"
#include "tidysq/util/calculate_length.h"
#include "tidysq/ProtoSequenceInputInterpreter.h"

namespace tidysq {
    template<typename INTERNAL, typename PROTO, bool SIMPLE>
    class ProtoSequenceInputInterpreter;

    template<typename INTERNAL, typename PROTO>
    class ProtoSequence {
        typename util::TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentStorageType content_;
    public:
        typedef typename PROTO::ProtoSequenceElementType                                    ElementType;
        typedef typename util::TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentStorageType    ContentStorageType;
        typedef typename util::TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentAccessType            AccessType;
        typedef typename util::TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentConstAccessType       ConstAccessType;

        explicit ProtoSequence(const ContentStorageType &content) :
                content_(content) {};

        explicit ProtoSequence(const LenSq length) :
                ProtoSequence(ContentStorageType(length)) {};

        ProtoSequence() :
                ProtoSequence(0) {};

        ProtoSequence(const std::initializer_list<ElementType> &list) :
                content_(ContentStorageType(list)) {};

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        };

        inline ConstAccessType operator[](const LenSq index) const {
            return content_[index];
        }

        [[nodiscard]] inline LenSq size() const {
            return content_.size();
        }

        [[nodiscard]] inline const ContentStorageType &content() const {
            return content_;
        }

        inline bool operator==(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return content_ == other.content_;
        }

        inline bool operator!=(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return !operator==(other);
        }

        template<bool ENABLED = std::is_same_v<PROTO, STRING_PT>>
        inline ProtoSequence<INTERNAL, PROTO>& operator+=(std::enable_if_t<ENABLED, const Letter &> letter) {
            content_ += letter;
            return *this;
        }

        template<bool ENABLED = std::is_same_v<PROTO, STRING_PT>>
        inline ProtoSequence<INTERNAL, PROTO>& operator+=(std::enable_if_t<!ENABLED, const Letter &> letter) {}

        template<bool ENABLED = std::is_same_v<PROTO, STRING_PT>>
        inline ProtoSequence<INTERNAL, PROTO>& operator+=(std::enable_if_t<ENABLED, const SimpleLetter &> letter) {
            content_ += letter;
            return *this;
        }

        template<bool ENABLED = std::is_same_v<PROTO, STRING_PT>>
        inline ProtoSequence<INTERNAL, PROTO>& operator+=(std::enable_if_t<!ENABLED, const SimpleLetter &> letter) {}

        template<bool SIMPLE>
        inline ProtoSequenceInputInterpreter<INTERNAL, PROTO, SIMPLE> content_interpreter(const Alphabet &alphabet) const {
            return ProtoSequenceInputInterpreter<INTERNAL, PROTO, SIMPLE>(content_.cbegin(), content_.cend(), alphabet);
        }
    };

    template<>
    inline ProtoSequence<STD_IT, STRING_PT>::ProtoSequence(const LenSq length) :
            ProtoSequence(ContentStorageType(length, ' ')) {}

    template<>
    inline ProtoSequence<RCPP_IT, STRING_PT>::ProtoSequence(const LenSq length) :
            ProtoSequence(ContentStorageType(length, ' ')) {}

    template<>
    inline bool ProtoSequence<RCPP_IT, RAWS_PT>::operator==(const ProtoSequence<RCPP_IT, RAWS_PT> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }

    template<>
    inline bool ProtoSequence<RCPP_IT, INTS_PT>::operator==(const ProtoSequence<RCPP_IT, INTS_PT> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }

    template<>
    inline bool ProtoSequence<RCPP_IT, STRINGS_PT>::operator==(const ProtoSequence<RCPP_IT, STRINGS_PT> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }
}

