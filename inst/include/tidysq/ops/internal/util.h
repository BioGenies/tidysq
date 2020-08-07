#ifndef TIDYSQ_UTIL_H
#define TIDYSQ_UTIL_H

#include <cmath>
#include "../../types/SequenceSTD.h"
#include "../../types/SequenceRCPP.h"


namespace tidysq::internal {
    template<InternalType INTERNAL>
    inline sizealph getAlphabetSize(const Alphabet<INTERNAL> &alphabet) {
        return ceil(log2((double) alphabet.size() + 2));
    }

    template<InternalType INTERNAL,
            ProtoType PROTO>
    inline lensq getPackedLength(const SequenceProto<INTERNAL, PROTO> &unpacked, const Alphabet<INTERNAL> &alphabet) {
        return (getAlphabetSize(alphabet) * unpacked.size() + 7) / 8;
    }

    template<InternalType INTERNAL>
    inline lensq getOriginalLength(const Sequence<INTERNAL> &packed);

    template<>
    inline lensq getOriginalLength<STD>(const Sequence<STD> &packed) {
        return packed.originalLength();
    }

    template<>
    inline lensq getOriginalLength<RCPP>(const Sequence<RCPP> &packed) {
        return packed.attr("original_length");
    }

    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserveSpaceForPacked(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        return Sequence<INTERNAL_OUT>(getPackedLength(unpacked, alphabet), unpacked.size());
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT,
            ProtoType PROTO_OUT>
    inline SequenceProto<INTERNAL_OUT, PROTO_OUT> reserveSpaceForUnpacked(const Sequence<INTERNAL_IN> &packed) {
        return SequenceProto<INTERNAL_OUT, PROTO_OUT>(getOriginalLength(packed));
    }

    template<InternalType INTERNAL>
    inline letvalue getNAValue(const Alphabet<INTERNAL> &alphabet) {
        return pow(2, getAlphabetSize(alphabet)) - 1;
    }

    template<InternalType INTERNAL>
    inline letvalue match(const std::string &letter, const Alphabet<INTERNAL> &alphabet, const letvalue &value) {
        for (letvalue i = 0; i < alphabet.size(); i++) {
            if (letter == alphabet[i]) {
                return i + 1;
            }
        }
        return value;
    }
}

#endif //TIDYSQ_UTIL_H
