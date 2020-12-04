#pragma once 

#include "tidysq/ops/Operation.h"
#include "tidysq/internal/pack.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationPack :
                public OperationVectorToVector<ProtoSq<INTERNAL_IN, PROTO_IN>, ProtoSequence<INTERNAL_IN, PROTO_IN>,
                                               Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>> {
            const Alphabet &alphabet_;
        public:
            explicit OperationPack(const Alphabet &alphabet) :
                    alphabet_(alphabet) {};

            inline Sq<INTERNAL_OUT> initialize_vector_out(const ProtoSq<INTERNAL_IN, PROTO_IN> &proto_sq, const LenSq from, const LenSq to) override {
                return Sq<INTERNAL_OUT>(to - from, alphabet_);
            }

            inline Sequence<INTERNAL_OUT> initialize_element_out(const ProtoSequence<INTERNAL_IN, PROTO_IN> &proto_sequence) override {
                return util::reserve_space_for_packed<INTERNAL_OUT>(proto_sequence.length(), alphabet_.alphabet_size());
            }

            inline void operator() (const ProtoSequence<INTERNAL_IN, PROTO_IN> &proto_sequence,
                                    Sequence<INTERNAL_OUT> &sequence) override {
                if (alphabet_.is_simple()) {
                    internal::pack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(proto_sequence, sequence, alphabet_);
                } else {
                    internal::pack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, false>(proto_sequence, sequence, alphabet_);
                }
            }
        };
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT = INTERNAL_IN>
    inline Sq<INTERNAL_OUT> pack(const ProtoSq<INTERNAL_IN, PROTO_IN> &proto_sq, const LenSq from, const LenSq to) {
        return sqapply(proto_sq, ops::OperationPack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(proto_sq.alphabet()), from, to);
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT = INTERNAL_IN>
    inline Sq<INTERNAL_OUT> pack(const ProtoSq<INTERNAL_IN, PROTO_IN> &proto_sq) {
        return pack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(proto_sq, 0, proto_sq.length());
    }


    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT = INTERNAL_IN>
    inline Sequence<INTERNAL_OUT> pack(const ProtoSequence<INTERNAL_IN, PROTO_IN> &proto_sequence, const Alphabet &alphabet) {
        return ops::OperationPack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(alphabet).
                template OperationVectorToVector<ProtoSq<INTERNAL_IN, PROTO_IN>, ProtoSequence<INTERNAL_IN, PROTO_IN>,
                                                      Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>>::operator() (proto_sequence);
    }
}