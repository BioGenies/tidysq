#pragma once 

#include "tidysq/ProtoSequence.h"

namespace tidysq::internal {
    inline Letter wrap_to_letter(const Letter &letter) {
        return letter;
    }

    inline Letter wrap_to_letter(const SimpleLetter &letter) {
        return {letter};
    }

    template<typename INTERNAL, bool SIMPLE>
    std::set<Letter> obtain_alphabet(const typename TypeBinder<INTERNAL, STRING_PT>::ProtoSqListConstructorType &x,
                                     const LenSq sample_size,
                                     const Letter &NA_letter,
                                     const bool ignore_case) {
        LenSq already_read = 0;
        std::set<Letter> letters = {};

        auto iter = x.begin();

        while (already_read < sample_size && iter != x.end()) {
            ProtoSequence<INTERNAL, STRING_PT> sequence((std::string(*iter)));
            auto interpreter = sequence.template content_interpreter<SIMPLE>(Alphabet(std::vector<Letter>{}, NA_letter));
            while (already_read < sample_size && !interpreter.reached_end()) {
                Letter letter = wrap_to_letter(interpreter.get_next_element());
                if (ignore_case && !isupper(letter[0])) {
                    letter[0] = std::toupper(letter[0]);
                }
                letters.insert(letter);
                already_read += 1;
            }
            iter++;
        }

        letters.erase(NA_letter);
        return letters;
    }
}