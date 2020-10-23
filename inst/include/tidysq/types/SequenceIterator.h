#ifndef TIDYSQ_SEQUENCEITERATOR_H
#define TIDYSQ_SEQUENCEITERATOR_H

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"
// #include "tidysq/ops/internal/util.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class SequenceIterator : public std::iterator<std::input_iterator_tag, ElementPacked> {
        typedef typename InternalTypeMapper<INTERNAL>::SequenceContentType ContentType;
        
        const ContentType& sequence_;
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
            
            inline ElementPacked operator*() const {
                ElementPacked ret = 0xff >> (8 - alph_size_);
                
                LenSq seq_lowest_bit_index = alph_size_ * pointer_;
                LenSq seq_highest_bit_index = seq_lowest_bit_index + alph_size_ - 1;
                LenSq seq_lowest_byte_index = seq_lowest_bit_index / 8;
                LenSq seq_highest_byte_index = seq_highest_bit_index / 8;
                unsigned short seq_lowest_bit_in_byte_index = seq_lowest_bit_index % 8;
                
                ret = ret &
                  ((sequence_[seq_lowest_byte_index] >> seq_lowest_bit_in_byte_index) |
                  (sequence_[seq_highest_byte_index] << (8 - seq_lowest_bit_in_byte_index)));
                
                return ret;
                // throw std::out_of_range("Write some error message, you fool!");
            }
            
            inline bool operator==(const SequenceIterator& other) const {
                return pointer_ == other.pointer_;
            }
          
            inline bool operator!=(const SequenceIterator& other) const {
                return !operator==(other);
            }
        
        friend class Sequence<INTERNAL>;
    };
}

#endif //TIDYSQ_SEQUENCEITERATOR_H
