#pragma once

#include "tidysq/types/general.h"

#include "tidysq/types/TypeMapper.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/types/SequenceIterator.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sequence {
        typename InternalTypeMapper<INTERNAL>::SequenceContentType content_;
        LenSq original_length_;
    public:
        typedef typename InternalTypeMapper<INTERNAL>::SequenceContentType ContentType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceElementType ElementType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceAccessType AccessType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceConstAccessType ConstAccessType;

        Sequence(const ContentType &content, const LenSq originalLength) :
                content_(content),
                original_length_(originalLength) {};

        Sequence(const LenSq contentLength, const LenSq originalLength) :
                Sequence(ContentType(contentLength), originalLength) {};

        Sequence() :
                Sequence(0, 0) {};

        Sequence(const Sequence &other) = default;

        Sequence(Sequence &&other) = default;

        Sequence& operator=(const Sequence &other) noexcept = default;

        Sequence& operator=(Sequence &&other) noexcept = default;

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator[](const LenSq index) const {
             return content_[index];
        }
        
        SequenceIterator<INTERNAL> begin(const Alphabet& alph) const {
            return SequenceIterator<INTERNAL>(content_, original_length_, alph);
        }
        
        SequenceIterator<INTERNAL> end(const Alphabet& alph) const {
            return SequenceIterator<INTERNAL>(content_, original_length_, alph, original_length_);
        }

        [[nodiscard]] inline LenSq originalLength() const {
            return original_length_;
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline ContentType content() const {
            return content_;
        }

        [[nodiscard]] inline bool operator==(const Sequence<INTERNAL> &other) const {
            return content_ == other.content_;
        }

        [[nodiscard]] inline bool operator!=(const Sequence<INTERNAL> &other) const {
            return !operator==(other);
        }

        void trim(const LenSq packed_length, const Alphabet &alphabet) {
            content_.erase(content_.begin() + internal::calculate_packed_internal_length(packed_length, alphabet), content_.end());
            original_length_ = packed_length;
        }

    };

    template<>
    inline bool Sequence<RCPP>::operator==(const Sequence<RCPP> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }
}
