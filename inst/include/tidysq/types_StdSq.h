#ifndef TIDYSQ_TYPES_STDSQ_H
#define TIDYSQ_TYPES_STDSQ_H

#include <vector>

#include "types_StdSequence.h"
#include "impl_Operation.h"
#include "sqapply.h"

namespace tidysq {
    class StdSq {
        std::vector<StdSequence> content_;
        Alphabet alphabet_;
    public:
        typedef StdSequence SequenceType;
        typedef Alphabet AlphabetType;

        StdSq(lensq length, AlphabetType alphabet) :
                content_(std::vector<StdSequence>(length)),
                alphabet_(alphabet) {};

        StdSequence &operator[] (lensq index) {
            return content_[index];
        }

        const StdSequence &operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] AlphabetType alphabet() const {
            return alphabet_;
        }

        template<typename TYPE_OUT>
        TYPE_OUT unpack() const {
            return sqapply<StdSq, TYPE_OUT, AlphabetType>(*this, ops::OperationUnpack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }
    };
}

#endif //TIDYSQ_TYPES_STDSQ_H
