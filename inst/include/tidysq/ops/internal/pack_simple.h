#ifndef TIDYSQ_PACK_SIMPLE_H
#define TIDYSQ_PACK_SIMPLE_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack2(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (*  interpreter      ) |
                                (*++interpreter << 2u) |
                                (*++interpreter << 4u) |
                                (*++interpreter << 6u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (*++interpreter      ) |
                                (*++interpreter << 2u) |
                                (*++interpreter << 4u) |
                                (*++interpreter << 6u);

            ++out_byte;
            ++interpreter;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack3(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (*  interpreter      ) |
                                (*++interpreter << 3u) |
                                ((tmp = *++interpreter) << 6u);

            if (++out_byte == packed.length()) return;
            
            packed[out_byte] =  (tmp            >> 2u) |
                                (*++interpreter << 1u) |
                                (*++interpreter << 4u) |
                                ((tmp = *++interpreter) << 7u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (tmp            >> 1u) |
                                (*++interpreter << 2u) |
                                (*++interpreter << 5u);
            
            ++interpreter;
            ++out_byte;
        }
    }


    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack4(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (*  interpreter      ) |
                                (*++interpreter << 4u);

            if (++out_byte == packed.length()) return;

            
            packed[out_byte] =  (*++interpreter      ) |
                                (*++interpreter << 4u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (*++interpreter      ) |
                                (*++interpreter << 4u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (*++interpreter      ) |
                                (*++interpreter << 4u);
            
            ++interpreter;
            ++out_byte;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT, bool SIMPLE>
    void pack5(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
               Sequence<INTERNAL_OUT> &packed,
               const Alphabet &alphabet) {
        LenSq out_byte = 0;
        auto interpreter = unpacked.template content_interpreter<SIMPLE>(alphabet);
        LetterValue tmp;
        while (!interpreter.reached_end()) {
            packed[out_byte] =  (*  interpreter       ) |
                                ((tmp = *++interpreter) << 5u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (tmp            >> 3u) |
                                (*++interpreter << 2u) |
                                ((tmp = *++interpreter) << 7u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (tmp            >> 1u) |
                                ((tmp = *++interpreter) << 4u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (tmp            >> 4u) |
                                (*++interpreter << 1u) |
                                ((tmp = *++interpreter) << 6u);

            if (++out_byte == packed.length()) return;

            packed[out_byte] =  (tmp            >> 2u) |
                                (*++interpreter << 3u);
            out_byte += 5;
        }
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    void pack_simple(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                     Sequence<INTERNAL_OUT> &packed,
                     const Alphabet &alphabet) {
        switch (alphabet.alphabet_size()) {
            case 2:
                pack2<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(unpacked, packed, alphabet);
                return;
            case 3:
                pack3<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(unpacked, packed, alphabet);
                return;
            case 4:
                pack4<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(unpacked, packed, alphabet);
                return;
            case 5:
                pack5<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(unpacked, packed, alphabet);
                return;
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACK_SIMPLE_H
