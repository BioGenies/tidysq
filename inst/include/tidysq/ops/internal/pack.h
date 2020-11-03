#ifndef TIDYSQ_PACK_H
#define TIDYSQ_PACK_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack2(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (interpreter.get_next()      ) |
                                (interpreter.get_next() << 2u) |
                                (interpreter.get_next() << 4u) |
                                (interpreter.get_next() << 6u);
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack3(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (interpreter.get_next()      ) |
                                (interpreter.get_next() << 3u) |
                                ((tmp = interpreter.get_next()) << 6u);

            if (++out_byte == packed.length()) break;
            
            packed[out_byte] =  (tmp                    >> 2u) |
                                (interpreter.get_next() << 1u) |
                                (interpreter.get_next() << 4u) |
                                ((tmp = interpreter.get_next()) << 7u);

            if (++out_byte == packed.length()) break;

            packed[out_byte] =  (tmp                    >> 1u) |
                                (interpreter.get_next() << 2u) |
                                (interpreter.get_next() << 5u);

            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }


    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack4(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (interpreter.get_next()      ) |
                                (interpreter.get_next() << 4u);
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack5(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (interpreter.get_next()       ) |
                                ((tmp = interpreter.get_next()) << 5u);

            if (++out_byte == packed.length()) break;

            packed[out_byte] =  (tmp                    >> 3u) |
                                (interpreter.get_next() << 2u) |
                                ((tmp = interpreter.get_next()) << 7u);

            if (++out_byte == packed.length()) break;

            packed[out_byte] =  (tmp                    >> 1u) |
                                ((tmp = interpreter.get_next()) << 4u);

            if (++out_byte == packed.length()) break;

            packed[out_byte] =  (tmp                    >> 4u) |
                                (interpreter.get_next() << 1u) |
                                ((tmp = interpreter.get_next()) << 6u);

            if (++out_byte == packed.length()) break;

            packed[out_byte] =  (tmp                    >> 2u) |
                                (interpreter.get_next() << 3u);
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
              Sequence<INTERNAL_OUT> &packed,
              const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                pack2<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, SIMPLE>(unpacked, packed, alphabet);
                return;
            case 3:
                pack3<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, SIMPLE>(unpacked, packed, alphabet);
                return;
            case 4:
                pack4<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, SIMPLE>(unpacked, packed, alphabet);
                return;
            case 5:
                pack5<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, SIMPLE>(unpacked, packed, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACK_H