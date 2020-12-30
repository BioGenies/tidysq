#pragma once

#include "tidysq/Alphabet.h"

namespace tidysq {
    template<typename INTERNAL, typename PROTO>
    class ProtoSequence;

    template<typename INTERNAL>
    class Sequence;
}


namespace tidysq::util {
    inline LenSq calculate_packed_internal_length(const LenSq unpacked_length, const AlphSize &alph_size) {
        return (alph_size * unpacked_length + 7) / 8;
    }

    template<typename INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserve_space_for_packed(const LenSq &unpacked_length,
                                                           const AlphSize &alph_size) {
        return Sequence<INTERNAL_OUT>(calculate_packed_internal_length(unpacked_length, alph_size), unpacked_length);
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> reserve_space_for_unpacked(const Sequence<INTERNAL_IN> &packed) {
        return ProtoSequence<INTERNAL_OUT, PROTO_OUT>(packed.original_length());
    }
}
