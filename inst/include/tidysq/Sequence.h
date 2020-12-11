#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/util/calculate_length.h"

namespace tidysq {
    template<typename INTERNAL>
    class Sequence {
        typename INTERNAL::SequenceContentStorageType content_;
        LenSq original_length_;
    public:
        typedef typename INTERNAL::SequenceContentStorageType   ContentStorageType;
        typedef typename INTERNAL::SequenceElementType          ElementType;
        typedef typename INTERNAL::SequenceContentAccessType           AccessType;
        typedef typename INTERNAL::SequenceContentConstAccessType      ConstAccessType;

    private:
        template<bool CONST>
        class GenericSequenceIterator : public std::iterator<std::bidirectional_iterator_tag, ElementPacked> {
        public:
            typedef std::conditional_t<CONST, const Sequence &, Sequence &> SequenceReference;
        protected:
            SequenceReference sequence_;
            const AlphSize alph_size_;
            LenSq pointer_;

        public:
            GenericSequenceIterator(SequenceReference sequence, const AlphSize &alph_size, const LenSq pointer);
            GenericSequenceIterator(SequenceReference sequence, const AlphSize &alph_size);
            GenericSequenceIterator(const GenericSequenceIterator<false> &other);

            GenericSequenceIterator& operator++();
            GenericSequenceIterator operator++(int);
            GenericSequenceIterator& operator--();
            GenericSequenceIterator operator--(int);

            // TODO: issue #26 (especially swap())
            ElementPacked operator*() const;
            bool operator==(const GenericSequenceIterator& other) const;
            bool operator!=(const GenericSequenceIterator& other) const;
            bool operator>(const GenericSequenceIterator& other) const;
            bool operator<(const GenericSequenceIterator& other) const;
            bool operator>=(const GenericSequenceIterator& other) const;
            bool operator<=(const GenericSequenceIterator& other) const;
            GenericSequenceIterator operator+(LenSq i) const;
            GenericSequenceIterator operator+(const GenericSequenceIterator& it) const;
            GenericSequenceIterator& operator+=(LenSq i);
            GenericSequenceIterator operator-(LenSq i) const;
            GenericSequenceIterator operator-(const GenericSequenceIterator& it) const;
            GenericSequenceIterator& operator-=(LenSq i);
            ElementPacked operator[](LenSq i);
            [[nodiscard]] LenSq index() const;

            template<bool ENABLED = !CONST>
            void assign(std::enable_if_t<ENABLED, const ElementPacked &> value);

            template<bool ENABLED = !CONST>
            void assign(std::enable_if_t<!ENABLED, const ElementPacked &> value);
        };

    public:
        typedef GenericSequenceIterator<true> const_iterator;
        typedef GenericSequenceIterator<false> iterator;

        Sequence(const ContentStorageType &content, const LenSq original_length) :
                content_(content),
                original_length_(original_length) {};

        Sequence(const LenSq content_length, const LenSq original_length) :
                Sequence(ContentStorageType(content_length), original_length) {};

        Sequence() :
                Sequence(0, 0) {};

        Sequence(const Sequence &other) = default;

        Sequence(Sequence &&other) noexcept = default;

        Sequence& operator=(const Sequence &other) = default;

        Sequence& operator=(Sequence &&other) noexcept = default;

        inline LetterValue operator[](const std::pair<LenSq, AlphSize> &index) const {
            ElementPacked ret = 0xffu >> (8u - std::get<1>(index));

            LenSq lowest_bit_index = std::get<1>(index) * std::get<0>(index);
            LenSq highest_bit_index = lowest_bit_index + std::get<1>(index) - 1;
            LenSq lowest_byte_index = lowest_bit_index / 8;
            LenSq highest_byte_index = highest_bit_index / 8;
            unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

            ret = ret &
                  ((content_[lowest_byte_index] >> lowest_bit_in_byte_index) |
                   (content_[highest_byte_index] << (8 - lowest_bit_in_byte_index)));

            return ret;
        }

        inline AccessType operator()(const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator()(const LenSq index) const {
            return content_[index];
        }
        
        iterator begin(const AlphSize& alph_size) {
            return iterator(*this, alph_size);
        }
        
        iterator end(const AlphSize& alph_size) {
            return iterator(*this, alph_size, original_length_);
        }

        const_iterator begin(const AlphSize& alph_size) const {
            return const_iterator(*this, alph_size);
        }

        const_iterator end(const AlphSize& alph_size) const {
            return const_iterator(*this, alph_size, original_length_);
        }

        const_iterator cbegin(const AlphSize& alph_size) const {
            return const_iterator(*this, alph_size);
        }

        const_iterator cend(const AlphSize& alph_size) const {
            return const_iterator(*this, alph_size, original_length_);
        }

        [[nodiscard]] inline LenSq original_length() const {
            return original_length_;
        }

        [[nodiscard]] inline LenSq size() const {
            return content_.size();
        }

        [[nodiscard]] inline ContentStorageType content() const {
            return content_;
        }

        [[nodiscard]] inline bool operator==(const Sequence<INTERNAL> &other) const {
            return content_ == other.content_;
        }

        [[nodiscard]] inline bool operator!=(const Sequence<INTERNAL> &other) const {
            return !operator==(other);
        }

        void trim(const LenSq packed_length, const Alphabet &alphabet) {
            content_.erase(content_.begin() + util::calculate_packed_internal_length(packed_length, alphabet.alphabet_size()), content_.end());
            original_length_ = packed_length;
        }
    };

    template<>
    inline bool Sequence<RCPP_IT>::operator==(const Sequence<RCPP_IT> &other) const {
        return Rcpp::is_true(Rcpp::all(content_ == other.content_));
    }

    //SEQUENCE<INTERNAL>::SEQUENCE_ITERATOR DEFINITION
    template<typename INTERNAL>
    template<bool CONST>
    inline Sequence<INTERNAL>::GenericSequenceIterator<CONST>::GenericSequenceIterator(
            SequenceReference sequence, const AlphSize &alph_size, const LenSq pointer) :
            sequence_(sequence),
            alph_size_(alph_size),
            pointer_(pointer) {}

    template<typename INTERNAL>
    template<bool CONST>
    inline Sequence<INTERNAL>::GenericSequenceIterator<CONST>::GenericSequenceIterator(
            SequenceReference sequence, const AlphSize &alph_size) :
            GenericSequenceIterator(sequence, alph_size, 0) {}


    template<typename INTERNAL>
    template<bool CONST>
    Sequence<INTERNAL>::GenericSequenceIterator<CONST>::GenericSequenceIterator(
            const Sequence::GenericSequenceIterator<false> &other) :
            sequence_(other.sequence_),
            alph_size_(other.alph_size_),
            pointer_(other.pointer_) {}

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST>
            &Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator++() {
        ++pointer_;
        return *this;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST>
            Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator++(int) {
        GenericSequenceIterator tmp(*this);
        operator++();
        return tmp;
    }
    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST>
            &Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator--() {
        --pointer_;
        return *this;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST>
            Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator--(int) {
        GenericSequenceIterator tmp(*this);
        operator--();
        return tmp;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline ElementPacked Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator*() const {
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

    template<typename INTERNAL>
    template<bool CONST>
    template<bool ENABLED>
    void Sequence<INTERNAL>::GenericSequenceIterator<CONST>::assign(std::enable_if_t<ENABLED, const ElementPacked &> value) {
        LenSq lowest_bit_index = alph_size_ * pointer_;
        LenSq highest_bit_index = lowest_bit_index + alph_size_ - 1;
        LenSq lowest_byte_index = lowest_bit_index / 8;
        LenSq highest_byte_index = highest_bit_index / 8;
        unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

        sequence_.content_[lowest_byte_index] = sequence_.content_[lowest_byte_index] |
                                                (value << lowest_bit_in_byte_index);
        if (highest_byte_index != lowest_byte_index) {
            sequence_.content_[highest_byte_index] = sequence_.content_[highest_byte_index] |
                                                     (value >> (8u - lowest_bit_in_byte_index));
        }
    }

    template<typename INTERNAL>
    template<bool CONST>
    template<bool ENABLED>
    void Sequence<INTERNAL>::GenericSequenceIterator<CONST>::assign(std::enable_if_t<!ENABLED, const ElementPacked &> value) {}

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator==(const GenericSequenceIterator &other) const {
        return pointer_ == other.pointer_;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator!=(const GenericSequenceIterator &other) const {
        return !operator==(other);
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator>(const GenericSequenceIterator &other) const {
        return pointer_ > other.pointer_;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator<(const GenericSequenceIterator &other) const {
        return pointer_ < other.pointer_;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator>=(const GenericSequenceIterator &other) const {
        return !operator<(other);
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline bool Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator<=(const GenericSequenceIterator &other) const {
        return !operator>(other);
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator+(LenSq i) const {
        GenericSequenceIterator tmp(*this);
        tmp += i;
        return tmp;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator+(const GenericSequenceIterator &it) const {
        return operator+(it.pointer_);
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> &Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator+=(LenSq i) {
        if (i + pointer_ > sequence_.original_length_)
            throw std::out_of_range("SequenceIterator tried to increment the pointer after its end.");
        pointer_ += i;
        return *this;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator-(LenSq i) const {
        GenericSequenceIterator tmp(*this);
        tmp -= i;
        return tmp;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator-(const GenericSequenceIterator &it) const {
        operator-(it.pointer_);
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline typename Sequence<INTERNAL>::template GenericSequenceIterator<CONST> &Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator-=(LenSq i) {
        if (i > pointer_) {
            throw std::out_of_range("SequenceIterator tried to decrement the pointer before its front.");
        }
        pointer_ -= i;
        return *this;
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline ElementPacked Sequence<INTERNAL>::GenericSequenceIterator<CONST>::operator[](LenSq i) {
        pointer_ += i;
        return operator*();
    }

    template<typename INTERNAL>
    template<bool CONST>
    inline LenSq Sequence<INTERNAL>::GenericSequenceIterator<CONST>::index() const {
        return pointer_;
    }
}
