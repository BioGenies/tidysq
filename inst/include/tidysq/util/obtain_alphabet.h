#pragma once

#include "tidysq/types/ProtoSequence.h"

namespace tidysq {
    namespace internal {
        inline Letter wrap_to_letter(const Letter &letter) {
            return letter;
        }

        inline Letter wrap_to_letter(const SimpleLetter &letter) {
            return {letter};
        }

        //TODO: currently there's no difference between both <I, S>::PSCT - fix it!
        template<InternalType INTERNAL, bool SIMPLE>
        std::set<Letter> obtain_alphabet(const typename TypeMapper<INTERNAL, STRING>::ProtoSqContentType &x,
                                         const LenSq sample_size,
                                         const Letter &NA_letter) {
            LenSq already_read = 0;
            std::set<Letter> letters = {};

            auto iter = x.begin();

            while (already_read < sample_size && iter != x.end()) {
                ProtoSequence<INTERNAL, STRING> sequence((char *) *iter);
                auto interpreter = sequence.template content_interpreter<SIMPLE>(Alphabet({}, NA_letter));
                while (already_read < sample_size && !interpreter.reached_end()) {
                    letters.insert(wrap_to_letter(interpreter.get_next_element()));
                    already_read += 1;
                }
                iter++;
            }

            letters.erase(NA_letter);

            return letters;
        }
    }

    namespace util {
        template<InternalType INTERNAL>
        Alphabet obtain_alphabet(const typename TypeMapper<INTERNAL, STRING>::ProtoSqContentType &x,
                                 const LenSq sample_size,
                                 const Letter &NA_letter) {

            std::set<Letter> letters;

            if (NA_letter.length() == 0) {
                throw std::invalid_argument("'NA_letter' should have at least one character!");
            } else if (NA_letter.length() == 1) {
                letters = tidysq::internal::obtain_alphabet<INTERNAL, true>(x, sample_size, NA_letter);
            } else {
                letters = tidysq::internal::obtain_alphabet<INTERNAL, false>(x, sample_size, NA_letter);
            }

            return Alphabet(convert_set_to_vector(letters),UNT, NA_letter);
        }
    }
}