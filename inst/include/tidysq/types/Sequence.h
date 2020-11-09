#pragma once

#include "tidysq/types/general.h"

#include "tidysq/types/TypeMapper.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq {
    template<InternalType INTERNAL>
    Sequence<INTERNAL> bite(
            const typename Sequence<INTERNAL>::ConstSequenceIterator &it, const std::vector<int> &indices);

    template<InternalType INTERNAL>
    class Sequence {
        typename InternalTypeMapper<INTERNAL>::SequenceContentType content_;
        LenSq original_length_;
    public:
        typedef typename InternalTypeMapper<INTERNAL>::SequenceContentType ContentType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceElementType ElementType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceAccessType AccessType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceConstAccessType ConstAccessType;

    private:
        class AbstractSequenceIterator : public std::iterator<std::bidirectional_iterator_tag, ElementPacked> {
        protected:
            Sequence &sequence_;
            const AlphSize alph_size_;
            LenSq pointer_;

        public:
            AbstractSequenceIterator(Sequence & sequence, const AlphSize &alph_size, const LenSq pointer);
            AbstractSequenceIterator(Sequence & sequence, const AlphSize &alph_size);

            AbstractSequenceIterator& operator++();
            AbstractSequenceIterator operator++(int);
            AbstractSequenceIterator& operator--();
            AbstractSequenceIterator operator--(int);

            // TODO: possibly implement swap()
            ElementPacked operator*() const;
            ElementPacked access(LenSq index);
            bool operator==(const AbstractSequenceIterator& other) const;
            bool operator!=(const AbstractSequenceIterator& other) const;
            bool operator>(const AbstractSequenceIterator& other) const;
            bool operator<(const AbstractSequenceIterator& other) const;
            bool operator>=(const AbstractSequenceIterator& other) const;
            bool operator<=(const AbstractSequenceIterator& other) const;
            AbstractSequenceIterator operator+(LenSq i) const;
            AbstractSequenceIterator operator+(const AbstractSequenceIterator& it) const;
            AbstractSequenceIterator& operator+=(LenSq i);
            AbstractSequenceIterator operator-(LenSq i) const;
            AbstractSequenceIterator operator-(const AbstractSequenceIterator& it) const;
            AbstractSequenceIterator& operator-=(LenSq i);
            ElementPacked operator[](LenSq i);
            [[nodiscard]] LenSq index() const;
        };

        class ConstSequenceIterator : public AbstractSequenceIterator {
        public:
            using AbstractSequenceIterator::AbstractSequenceIterator;

            friend Sequence<INTERNAL> bite<INTERNAL>(
                    const ConstSequenceIterator &it, const std::vector<int> &indices);
        };

        class SequenceIterator : public AbstractSequenceIterator {
        public:
            using AbstractSequenceIterator::AbstractSequenceIterator;

            void assign(const ElementPacked &value);
        };

    public:
        Sequence(const ContentType &content, const LenSq originalLength) :
                content_(content),
                original_length_(originalLength) {};

        Sequence(const LenSq contentLength, const LenSq originalLength) :
                Sequence(ContentType(contentLength), originalLength) {};

        Sequence() :
                Sequence(0, 0) {};

        Sequence(const Sequence &other) = default;

        Sequence(Sequence &&other) noexcept = default;

        Sequence& operator=(const Sequence &other) = default;

        Sequence& operator=(Sequence &&other) noexcept = default;

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator[](const LenSq index) const {
             return content_[index];
        }
        
        SequenceIterator begin(const AlphSize& alph_size) {
            return SequenceIterator(*this, alph_size);
        }
        
        SequenceIterator end(const AlphSize& alph_size) {
            return SequenceIterator(*this, alph_size, original_length_);
        }

        ConstSequenceIterator cbegin(const AlphSize& alph_size) const {
            return ConstSequenceIterator(*this, alph_size);
        }

        SequenceIterator cend(const AlphSize& alph_size) const {
            return SequenceIterator(*this, alph_size, original_length_);
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
            content_.erase(content_.begin() + internal::calculate_packed_internal_length(packed_length, alphabet.alphabet_size()), content_.end());
            original_length_ = packed_length;
        }
    };

    template<>
    inline bool Sequence<RCPP>::operator==(const Sequence<RCPP> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }

    //SEQUENCE<INTERNAL>::SEQUENCE_ITERATOR DEFINITION
    template<InternalType INTERNAL>
    inline Sequence<INTERNAL>::AbstractSequenceIterator::AbstractSequenceIterator(Sequence &sequence, const AlphSize &alph_size, const LenSq pointer) :
            sequence_(sequence),
            alph_size_(alph_size),
            pointer_(pointer) {}

    template<InternalType INTERNAL>
    inline Sequence<INTERNAL>::AbstractSequenceIterator::AbstractSequenceIterator(Sequence &sequence, const AlphSize &alph_size) :
            AbstractSequenceIterator(sequence, alph_size, 0) {}

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator &Sequence<INTERNAL>::AbstractSequenceIterator::operator++() {
        ++pointer_;
        return *this;
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator++(int) {
        AbstractSequenceIterator tmp(*this);
        operator++();
        return tmp;
    }
    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator &Sequence<INTERNAL>::AbstractSequenceIterator::operator--() {
        --pointer_;
        return *this;
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator--(int) {
        AbstractSequenceIterator tmp(*this);
        operator--();
        return tmp;
    }

    template<InternalType INTERNAL>
    inline ElementPacked Sequence<INTERNAL>::AbstractSequenceIterator::operator*() const {
        ElementPacked ret = 0xffu >> (8u - alph_size_);

        LenSq lowest_bit_index = alph_size_ * pointer_;
        LenSq highest_bit_index = lowest_bit_index + alph_size_ - 1;
        LenSq lowest_byte_index = lowest_bit_index / 8;
        LenSq highest_byte_index = highest_bit_index / 8;
        unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

        ret = ret &
              ((sequence_.content_[lowest_byte_index] >> lowest_bit_in_byte_index) |
               (sequence_.content_[highest_byte_index] << (8 - lowest_bit_in_byte_index)));

        return ret;
    }

    template<InternalType INTERNAL>
    inline ElementPacked Sequence<INTERNAL>::AbstractSequenceIterator::access(LenSq index) {
        pointer_ = index;
        return operator*();
    }

    template<InternalType INTERNAL>
    void Sequence<INTERNAL>::SequenceIterator::assign(const ElementPacked &value) {
        LenSq lowest_bit_index = AbstractSequenceIterator::alph_size_ * AbstractSequenceIterator::pointer_;
        LenSq highest_bit_index = lowest_bit_index + AbstractSequenceIterator::alph_size_ - 1;
        LenSq lowest_byte_index = lowest_bit_index / 8;
        LenSq highest_byte_index = highest_bit_index / 8;
        unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

        AbstractSequenceIterator::sequence_.content_[lowest_byte_index] = AbstractSequenceIterator::sequence_.content_[lowest_byte_index] |
                                                (value << lowest_bit_in_byte_index);
        if (highest_byte_index != lowest_byte_index) {
            AbstractSequenceIterator::sequence_.content_[highest_byte_index] = AbstractSequenceIterator::sequence_.content_[highest_byte_index] |
                                                     (value >> (8u - lowest_bit_in_byte_index));
        }
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator==(const AbstractSequenceIterator &other) const {
        return pointer_ == other.pointer_;
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator!=(const AbstractSequenceIterator &other) const {
        return !operator==(other);
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator>(const AbstractSequenceIterator &other) const {
        return pointer_ > other.pointer_;
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator<(const AbstractSequenceIterator &other) const {
        return pointer_ < other.pointer_;
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator>=(const AbstractSequenceIterator &other) const {
        return !operator<(other);
    }

    template<InternalType INTERNAL>
    inline bool Sequence<INTERNAL>::AbstractSequenceIterator::operator<=(const AbstractSequenceIterator &other) const {
        return !operator>(other);
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator+(LenSq i) const {
        AbstractSequenceIterator tmp(*this);
        tmp += i;
        return tmp;
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator+(const AbstractSequenceIterator &it) const {
        return operator+(it.pointer_);
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator &Sequence<INTERNAL>::AbstractSequenceIterator::operator+=(LenSq i) {
        if (i + pointer_ > sequence_.originalLength_)
            throw std::out_of_range("SequenceIterator tried to increment the pointer after its end.");
        pointer_ += i;
        return *this;
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator-(LenSq i) const {
        AbstractSequenceIterator tmp(*this);
        tmp -= i;
        return tmp;
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator Sequence<INTERNAL>::AbstractSequenceIterator::operator-(const AbstractSequenceIterator &it) const {
        operator-(it.pointer_);
    }

    template<InternalType INTERNAL>
    inline typename Sequence<INTERNAL>::AbstractSequenceIterator &Sequence<INTERNAL>::AbstractSequenceIterator::operator-=(LenSq i) {
        if (i > pointer_) {
            throw std::out_of_range("SequenceIterator tried to decrement the pointer before its front.");
        }
        pointer_ -= i;
        return *this;
    }

    template<InternalType INTERNAL>
    inline ElementPacked Sequence<INTERNAL>::AbstractSequenceIterator::operator[](LenSq i) {
        pointer_ += i;
        return operator*();
    }

    template<InternalType INTERNAL>
    inline LenSq Sequence<INTERNAL>::AbstractSequenceIterator::index() const {
        return pointer_;
    }
}
