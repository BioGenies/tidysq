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
    std::set<Letter> obtain_alphabet(const typename TypeMapper<INTERNAL, STRING_PT>::ProtoSqContentType &x,
                                     const LenSq sample_size,
                                     const Letter &NA_letter,
                                     const bool ignore_case) {
        LenSq already_read = 0;
        std::set<Letter> letters = {};

        auto iter = x.begin();

        while (already_read < sample_size && iter != x.end()) {
            ProtoSequence<INTERNAL, STRING_PT> sequence((char *) *iter);
            auto interpreter = sequence.template content_interpreter<SIMPLE>(Alphabet(std::vector<Letter>{}, NA_letter));
            while (already_read < sample_size && !interpreter.reached_end()) {
                letters.insert(wrap_to_letter(interpreter.get_next_element()));
                already_read += 1;
            }
            iter++;
        }

        letters.erase(NA_letter);
        if (ignore_case) {
            std::set<Letter> letters_to_erase = {};
            for (const auto &letter : letters) {
                if (isalpha(letter[0]) && !isupper(letter[0]))
                    letters_to_erase.insert(letter);
            }
            for (const auto &letter : letters_to_erase) {
                letters.erase(letter);
            }
        }
        return letters;
    }
}