#ifndef TIDYSQ_PACKSTRING_H
#define TIDYSQ_PACKSTRING_H

#include "packSTRING_STANDARD.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packSTRING(const SequenceProto<INTERNAL_IN, STRING> &unpacked,
                                      const Alphabet &alphabet) {
        if (alphabet.isSimple()) {
            return packSTRING_STANDARD<INTERNAL_OUT>(unpacked, alphabet);
        }
        return Sequence<INTERNAL_OUT>(0, 0); // TODO: implement it!
    }
}

#endif //TIDYSQ_PACKSTRING_H
