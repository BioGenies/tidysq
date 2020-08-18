#ifndef TIDYSQ_UNPACKSTRING_H
#define TIDYSQ_UNPACKSTRING_H

#include "../../types/general.h"
#include "../../types/SequenceProtoSTD.h"
#include "../../types/SequenceProtoRCPP.h"
#include "unpackSTRING_STANDARD.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN>
    SequenceProto<STD, STRING> unpackSTRING_STD(const Sequence<INTERNAL_IN> &packed,
                                                const Alphabet<INTERNAL_IN> &alphabet) {
        if (alphabet.isSimple()) {
            return unpackSTRING_STANDARD<INTERNAL_IN>(packed, Alphabet<STD>(alphabet));
        }
        return SequenceProto<STD, STRING>(); // TODO: implement it!
    }

    template<InternalType INTERNAL_IN>
    SequenceProto<RCPP, STRING> unpackSTRING_RCPP(const Sequence<INTERNAL_IN> &packed,
                                                 const Alphabet<INTERNAL_IN> &alphabet) {
        return SequenceProto<RCPP, STRING>::steal(unpackSTRING_STD(packed, alphabet));
    }
}


#endif //TIDYSQ_UNPACKSTRING_H
