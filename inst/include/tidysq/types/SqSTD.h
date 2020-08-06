#ifndef TIDYSQ_SQSTD_H
#define TIDYSQ_SQSTD_H

#include <vector>

#include "general.h"
#include "SequenceSTD.h"
#include "AlphabetSTD.h"
#include "../impl_Operation.h"
#include "../sqapply.h"

namespace tidysq {
    template<>
    class Sq<STD> {
        std::vector<Sequence<STD>> content_;
        Alphabet<STD> alphabet_;
    public:
        typedef Sequence<STD> SequenceType;
        typedef Alphabet<STD> AlphabetType;

        Sq(lensq length, AlphabetType alphabet) :
                content_(std::vector<Sequence<STD>>(length)),
                alphabet_(alphabet) {};

        SequenceType &operator[] (lensq index) {
            return content_[index];
        }

        const SequenceType &operator[] (lensq index) const {
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
            return sqapply<Sq<STD>, TYPE_OUT, AlphabetType>(*this, ops::OperationUnpack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }
    };
}

#endif //TIDYSQ_SQSTD_H
