#pragma once

#include <map>

#include "tidysq/Sq.h"
#include "tidysq/ops/pack.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> substitute_letters(const Sequence<INTERNAL> &sequence,
                                          const std::unordered_map<Letter, Letter> &encoding,
                                          const Alphabet &alph,
                                          const Alphabet &dest_alph) {
        // TODO: maybe we can do something to avoid unpacking (and not lose performance)
        ProtoSequence<INTERNAL, STRINGS_PT> unpacked =
                util::reserve_space_for_unpacked<INTERNAL, INTERNAL, STRINGS_PT>(sequence);
        internal::unpack_common(sequence, unpacked, alph);

        // We have content as vector of strings, so that it's easier to swap them with encoding map
        ProtoSequence<INTERNAL, STRINGS_PT> encoded(unpacked.content());
        for (LenSq index = 0; index < encoded.length(); ++index) {
            Letter letter = alph[sequence[{index, alph.alphabet_size()}]];
            if (encoding.count(letter) > 0) {
                encoded[index] = encoding.at(letter);
            }
        }

        return pack<INTERNAL, STRINGS_PT, INTERNAL>(encoded, dest_alph);
    }

    template<typename INTERNAL>
    Sq<INTERNAL> substitute_letters(const Sq<INTERNAL> &sq,
                                    const std::unordered_map<Letter, Letter> &encoding) {
        const Alphabet &alph = sq.alphabet();

        // TODO: move the chunk below somewhere into Alphabet class
        std::vector<Letter> dest_letters;
        // We cannot rely on range-based loop, because it isn't ordered alphabetical by key. Thus it's necessary
        // to do manual looping. It should work as long as Alphabet class doesn't change to much.
        // TODO: come up with an idea for ordered looping that doesn't rely on this exact implementation of Alphabet
        for (LetterValue key = 0; key < alph.length(); ++key) {
            const Letter letter = alph[key];
            if (encoding.count(letter) == 0) {
                if (std::none_of(dest_letters.begin(), dest_letters.end(), [=](const Letter &other) {
                    return letter == other;
                })) {
                    dest_letters.push_back(letter);
                }
            } else if (std::none_of(dest_letters.begin(), dest_letters.end(), [=](const Letter &other) {
                return encoding.at(letter) == other;
            })) {
                dest_letters.push_back(encoding.at(letter));
            }
        }
        Alphabet dest_alph(dest_letters, ATP, alph.NA_letter());

        // If encoding doesn't reduce the number of letters, we can get away with just changing alphabet
        // Else we have to repack whole Sq object with new alphabet
        Sq<INTERNAL> ret(sq.length(), dest_alph);

        // First we have to extract all values from inside map (sadly, there's no easy method like in Java)
        std::vector<Letter> values;
        std::set<Letter> unique_values;
        for (const auto& entry : encoding) {
            values.push_back(entry.second);
            unique_values.insert(entry.second);
        }

        const bool need_repacking = std::any_of(values.cbegin(), values.cend(), [=](const Letter &value) {
            return alph.contains(value);
        }) || values.size() != unique_values.size();

        for (LenSq i = 0; i < sq.length(); ++i) {
            if (!need_repacking) {
                ret[i] = sq[i].get();
            } else {
                ret[i] = substitute_letters<INTERNAL>(sq[i].get(), encoding, alph, dest_alph);
            }
        }
        return ret;
    }
}
