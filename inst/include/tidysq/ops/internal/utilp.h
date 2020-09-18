#ifndef TIDYSQ_UTILP_H
#define TIDYSQ_UTILP_H

#include "tidysq/types/Sequence.h"
#include "tidysq/types/SequenceProto.h"

namespace tidysq::internal {
    template<ProtoType PROTO_OUT>
    auto matchLetter(letvalue value, const Alphabet &alphabet) -> typename SequenceProto<STD, PROTO_OUT>::ElementType;

    template<>
    inline ElemRaws matchLetter<RAWS>(const letvalue value, const Alphabet &alphabet) {
        return value;
    }
}

#endif //TIDYSQ_UTILP_H
