#ifndef TIDYSQ_TYPES_SQ_H
#define TIDYSQ_TYPES_SQ_H

#include <vector>

#include "types_Sequence.h"
#include "impl_Operation.h"
#include "sqapply.h"

namespace tidysq {
    class Sq {
        std::vector<Sequence> content_;
        Alphabet alphabet_;
    public:
        typedef Sequence SequenceType;
        typedef Alphabet AlphabetType;

        Sq(lensq length, AlphabetType alphabet) :
                content_(std::vector<Sequence>(length)),
                alphabet_(alphabet) {};

        Sequence &operator[] (lensq index) {
            return content_[index];
        }

        const Sequence &operator[] (lensq index) const {
            return content_[index];
        }

        lensq length() const {
            return content_.size();
        }

        AlphabetType alphabet() const {
            return alphabet_;
        }

        template<typename TYPE_OUT>
        TYPE_OUT unpack() const {
            return sqapply<Sq, TYPE_OUT, AlphabetType>(*this, ops::OperationUnpack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }
    };
}

#endif //TIDYSQ_TYPES_SQ_H
