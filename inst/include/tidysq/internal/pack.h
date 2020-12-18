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
        while (!interpreter.reached_end()) {
            packed(out_byte) = (interpreter.get_next_value()      ) |
                               (interpreter.get_next_value() << 2u) |
                               (interpreter.get_next_value() << 4u) |
                               (interpreter.get_next_value() << 6u) ;
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
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed(out_byte)  = (interpreter.get_next_value()      ) |
                                (interpreter.get_next_value() << 3u) ;
                         tmp  =  interpreter.get_next_value()        ;
            packed(out_byte) |= (tmp                          << 6u) ;

            if (++out_byte == packed.size()) break;
            
            packed(out_byte)  = (tmp                          >> 2u) |
                                (interpreter.get_next_value() << 1u) |
                                (interpreter.get_next_value() << 4u) ;
                         tmp  = interpreter.get_next_value()         ;
            packed(out_byte) |= (tmp                          << 7u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                          >> 1u) |
                                (interpreter.get_next_value() << 2u) |
                                (interpreter.get_next_value() << 5u) ;

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
        while (!interpreter.reached_end()) {
            packed(out_byte) = (interpreter.get_next_value()      ) |
                               (interpreter.get_next_value() << 4u) ;
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
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed(out_byte)  = (interpreter.get_next_value()         ) ;
                         tmp  =  interpreter.get_next_value()           ;
            packed(out_byte) |= (tmp                             << 5u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                             >> 3u) |
                                (interpreter.get_next_value()    << 2u) ;
                         tmp  =  interpreter.get_next_value()           ;
            packed(out_byte) |= (tmp                             << 7u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                             >> 1u) ;
                         tmp  =  interpreter.get_next_value()           ;
            packed(out_byte) |= (tmp                             << 4u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                             >> 4u) |
                                (interpreter.get_next_value()    << 1u) ;
                         tmp  =  interpreter.get_next_value()           ;
            packed(out_byte) |= (tmp                             << 6u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                             >> 2u) |
                                (interpreter.get_next_value()    << 3u) ;
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
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed(out_byte)  = (interpreter.get_next_value()      ) ;
                         tmp  =  interpreter.get_next_value()        ;
            packed(out_byte) |= (tmp                          << 6u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                          >> 2u) ;
                         tmp  = interpreter.get_next_value()         ;
            packed(out_byte) |= (tmp                          << 4u) ;

            if (++out_byte == packed.size()) break;

            packed(out_byte)  = (tmp                          >> 4u) |
                                (interpreter.get_next_value() << 2u) ;
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
