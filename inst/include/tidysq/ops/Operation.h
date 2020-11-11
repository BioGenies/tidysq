#ifndef TIDYSQ_OPERATION_H
#define TIDYSQ_OPERATION_H

#include "tidysq/Alphabet.h"

namespace tidysq::ops {
    template<typename SEQUENCE_IN,
            typename SEQUENCE_OUT>
    class OperationSq {
    public:
        virtual SEQUENCE_OUT operator() (const SEQUENCE_IN &sequence, const Alphabet &alphabet) const = 0;
    };
}

#endif //TIDYSQ_OPERATION_H
