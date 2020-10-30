#ifndef TIDYSQ_SEQUENCEITERATOR_H
#define TIDYSQ_SEQUENCEITERATOR_H

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class SequenceIterator : public std::iterator<std::bidirectional_iterator_tag, ElementPacked> {
        typedef typename InternalTypeMapper<INTERNAL>::SequenceContentType ContentType;
        
        ContentType sequence_;
        const AlphSize alph_size_;
        const LenSq& originalLength_;
        LenSq pointer_;
        
        SequenceIterator(const ContentType& content, const LenSq& length, const Alphabet& alph) :
                sequence_{content}, alph_size_{alph.alphabet_size()}, originalLength_{length}, pointer_{0} {}
        SequenceIterator(const ContentType& content, const LenSq& length, const Alphabet& alph, LenSq pointer) :
                sequence_{content}, alph_size_{alph.alphabet_size()}, originalLength_{length}, pointer_{pointer} {}
        
    public:
        inline SequenceIterator& operator++() {
            ++pointer_;
            return *this;
        }

        inline SequenceIterator operator++(int) {
            SequenceIterator tmp(*this);
            operator++();
            return tmp;
        }
        inline SequenceIterator& operator--() {
            --pointer_;
            return *this;
        }

        inline SequenceIterator operator--(int) {
            SequenceIterator tmp(*this);
            operator--();
            return tmp;
        }

        // TODO: possibly implement swap()

        inline ElementPacked operator*() const {
            ElementPacked ret = 0xffu >> (8u - alph_size_);

            LenSq lowest_bit_index = alph_size_ * pointer_;
            LenSq highest_bit_index = lowest_bit_index + alph_size_ - 1;
            LenSq lowest_byte_index = lowest_bit_index / 8;
            LenSq highest_byte_index = highest_bit_index / 8;
            unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

            ret = ret &
              ((sequence_[lowest_byte_index] >> lowest_bit_in_byte_index) |
               (sequence_[highest_byte_index] << (8 - lowest_bit_in_byte_index)));

            return ret;
        }

        inline ElementPacked access(const LenSq index) {
            pointer_ = index;
            return operator*();
        }

        inline void assign(const ElementPacked &value) {
            LenSq lowest_bit_index = alph_size_ * pointer_;
            LenSq highest_bit_index = lowest_bit_index + alph_size_ - 1;
            LenSq lowest_byte_index = lowest_bit_index / 8;
            LenSq highest_byte_index = highest_bit_index / 8;
            unsigned short lowest_bit_in_byte_index = lowest_bit_index % 8;

            sequence_[lowest_byte_index] = sequence_[lowest_byte_index] |
                    (value << lowest_bit_in_byte_index);
            if (highest_byte_index != lowest_byte_index) {
                sequence_[highest_byte_index] = sequence_[highest_byte_index] |
                        (value >> (8u - lowest_bit_in_byte_index));
            }
        }
            
        inline bool operator==(const SequenceIterator& other) const {
            return pointer_ == other.pointer_;
        }

        inline bool operator!=(const SequenceIterator& other) const {
            return !operator==(other);
        }

        inline bool operator>(const SequenceIterator& other) const {
            return pointer_ > other.pointer_;
        }

        inline bool operator<(const SequenceIterator& other) const {
            return pointer_ < other.pointer_;
        }

        inline bool operator>=(const SequenceIterator& other) const {
            return !operator<(other);
        }

        inline bool operator<=(const SequenceIterator& other) const {
            return !operator>(other);
        }

        inline SequenceIterator operator+(LenSq i) const {
            SequenceIterator tmp(*this);
            tmp += i;
            return tmp;
        }

        inline SequenceIterator operator+(const SequenceIterator& it) const {
            return operator+(it.pointer_);
        }

        inline SequenceIterator& operator+=(LenSq i) {
            if (i + pointer_ > originalLength_)
                throw std::out_of_range("SequenceIterator tried to increment the pointer after its end.");
            pointer_ += i;
            return *this;
        }

        inline SequenceIterator operator-(LenSq i) const {
            SequenceIterator tmp(*this);
            tmp -= i;
            return tmp;
        }

        inline SequenceIterator operator-(const SequenceIterator& it) const {
            operator-(it.pointer_);
        }

        inline SequenceIterator& operator-=(LenSq i) {
            if (i > pointer_) {
                throw std::out_of_range("SequenceIterator tried to decrement the pointer before its front.");
            }
            pointer_ -= i;
            return *this;
        }

        inline ElementPacked operator[](LenSq i) {
            pointer_ += i;
            return operator*();
        }

        friend class Sequence<INTERNAL>;
    };
}

#endif //TIDYSQ_SEQUENCEITERATOR_H
