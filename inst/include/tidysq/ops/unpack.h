#pragma once

#include "tidysq/internal/unpack_common.h"
#include "tidysq/internal/unpack_multichar.h"
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
        class OperationUnpack :
                public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                ProtoSq<INTERNAL_OUT, PROTO_OUT>, ProtoSequence<INTERNAL_OUT, PROTO_OUT>> {
            const Alphabet &alphabet_;
        public:
            explicit OperationUnpack(const Alphabet &alphabet) :
                    alphabet_(alphabet) {};

            inline ProtoSq<INTERNAL_OUT, PROTO_OUT> initialize_vector_out(const Sq<INTERNAL_IN> &sq, const LenSq from, const LenSq to) override {
                return ProtoSq<INTERNAL_OUT, PROTO_OUT>(to - from, alphabet_);
            }

            inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return util::reserve_space_for_unpacked<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(sequence);
            }

            inline void operator()(const Sequence<INTERNAL_IN> &sequence,
                       ProtoSequence<INTERNAL_OUT, PROTO_OUT> &proto_sequence) override {
                internal::unpack_common<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(sequence, proto_sequence, alphabet_);
            }
        };

        template<typename INTERNAL_IN, typename INTERNAL_OUT>
        class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, STRING_PT> :
                public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                ProtoSq<INTERNAL_OUT, STRING_PT>, ProtoSequence<INTERNAL_OUT, STRING_PT>> {
            const Alphabet &alphabet_;
            public:
            explicit OperationUnpack(const Alphabet &alphabet) :
                    alphabet_(alphabet) {};

            inline ProtoSq<INTERNAL_OUT, STRING_PT> initialize_vector_out(const Sq<INTERNAL_IN> &sq, const LenSq from, const LenSq to) override {
                return ProtoSq<INTERNAL_OUT, STRING_PT>(to - from, alphabet_);
            }

            ProtoSequence<INTERNAL_OUT, STRING_PT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                if (alphabet_.is_simple()) {
                    return util::reserve_space_for_unpacked<INTERNAL_IN, INTERNAL_OUT, STRING_PT>(sequence);
                } else {
                    return {};
                }
            }

            inline void operator()(const Sequence<INTERNAL_IN> &sequence,
                       ProtoSequence<INTERNAL_OUT, STRING_PT> &proto_sequence) override {
                if (alphabet_.is_simple()) {
                    internal::unpack_common<INTERNAL_IN, INTERNAL_OUT, STRING_PT>(sequence, proto_sequence, alphabet_);
                } else {
                    internal::unpack_multichar_string<INTERNAL_IN, INTERNAL_OUT>(sequence, proto_sequence, alphabet_);
                }
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    inline ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack(const Sq<INTERNAL_IN> &sq, const LenSq from, const LenSq to) {
        return sqapply(sq, ops::OperationUnpack<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(sq.alphabet()), from, to);
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    inline ProtoSq<INTERNAL_OUT, PROTO_OUT> unpack(const Sq<INTERNAL_IN> &sq) {
        return unpack<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(sq, 0, sq.size());
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT, typename PROTO_OUT>
    inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> unpack(const Sequence<INTERNAL_IN> &sequence, const Alphabet &alphabet) {
        return ops::OperationUnpack<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(alphabet).
                template OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                ProtoSq<INTERNAL_OUT, PROTO_OUT>, ProtoSequence<INTERNAL_OUT, PROTO_OUT>>::operator()(sequence);
    }
}