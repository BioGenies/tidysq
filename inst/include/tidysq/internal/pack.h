#pragma once

#include <stdexcept>

#include "tidysq/Alphabet.h"
#include "tidysq/Sequence.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/util/calculate_length.h"

namespace tidysq::internal {
    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
    void pack2(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue v1, v2, v3, v4;
        while (!interpreter.reached_end()) {
            v1 = interpreter.get_next_value();
            v2 = interpreter.get_next_value() << 2u;
            v3 = interpreter.get_next_value() << 4u;
            v4 = interpreter.get_next_value() << 6u;
            packed(out_byte) = v1 | v2 | v3 | v4 ;
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
    void pack3(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue v1, v2, v3, v4;
        while (!interpreter.reached_end()) {
            v1 = interpreter.get_next_value();
            v2 = interpreter.get_next_value() << 3u;
            v3 = interpreter.get_next_value();
            packed(out_byte)  = v1 | v2 | (v3 << 6u) ;
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value() << 1u;
            v2 = interpreter.get_next_value() << 4u;
            v4 = interpreter.get_next_value();
            packed(out_byte)  = (v3 >> 2u) | v1 | v2 | (v4 << 7u) ;
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value() << 2u;
            v2 = interpreter.get_next_value() << 5u;
            packed(out_byte)  = (v4 >> 1u) | v1 | v2 ;
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }


    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
    void pack4(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue v1, v2;
        while (!interpreter.reached_end()) {
            v1 = interpreter.get_next_value();
            v2 = interpreter.get_next_value() << 4u;
            packed(out_byte) = v1 | v2 ;
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
    void pack5(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue v1, v2, v3;
        while (!interpreter.reached_end()) {
            v1 = interpreter.get_next_value();
            v2 = interpreter.get_next_value();
            packed(out_byte)  = v1 | (v2 << 5u) ;
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value() << 2u;
            v3 = interpreter.get_next_value();
            packed(out_byte)  = (v2 >> 3u) | v1 | (v3 << 7u) ;
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value();
            packed(out_byte)  = (v3 >> 1u) | (v1 << 4u) ;
            if (++out_byte == packed.size()) break;

            v2 = interpreter.get_next_value() << 1u;
            v3 = interpreter.get_next_value();
            packed(out_byte)  = (v1 >> 4u) | v2 | (v3 << 6u) ;
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value() << 3u;
            packed(out_byte)  = (v3 >> 2u) | v1 ;
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
    void pack6(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue v1, v2;
        while (!interpreter.reached_end()) {
            v1 = interpreter.get_next_value();
            v2 = interpreter.get_next_value();
            packed(out_byte)  = v1 | (v2 << 6u);
            if (++out_byte == packed.size()) break;

            v1 = interpreter.get_next_value();
            packed(out_byte)  = (v2 >> 2u) | (v1 << 4u) ;
            if (++out_byte == packed.size()) break;

            v2 = interpreter.get_next_value() << 2u;
            packed(out_byte)  = (v1 >> 4u) | v2 ;
            ++out_byte;
        }
        packed.trim(interpreter.interpreted_letters(), alphabet);
    }

    template<typename INTERNAL_IN, typename PROTO_IN, typename INTERNAL_OUT, bool SIMPLE>
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
            case 6:
                pack6<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, SIMPLE>(unpacked, packed, alphabet);
                return;
            default:
                throw std::invalid_argument(std::string("\"alphabet\" has invalid alphabet size - it is ") +
                                              std::to_string(alphabet.alphabet_size()) +
                                              " but it should be between 2 and 6 inclusive");
        }
    }
}
