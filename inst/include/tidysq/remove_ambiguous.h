#pragma once

#include "tidysq/exports.h"

namespace tidysq {
    template<InternalType INTERNAL>
    ProtoSequence<INTERNAL, STRING> remove_ambiguous(const Sequence<INTERNAL> &sequence,
                                                     const Alphabet &alph,
                                                     const Alphabet &dest_alph,
                                                     const bool by_letter) {
        // TODO: extract functions out of these ifs, maybe?
        if (by_letter) {
            std::string selected_letters;
            for (auto it = sequence.cbegin(alph.alphabet_size()); it != sequence.cend(alph.alphabet_size()); ++it) {
                // We can suppose that Letters are simple, because only AMI_EXT, DNA_EXT and RNA_EXT are valid Sq objects
                Letter letter = alph[*it];
                if (dest_alph.contains(letter)) {
                    selected_letters += letter;
                }
            }
            return ProtoSequence<INTERNAL, STRING>{selected_letters};
        } else {
            // We check if all letters in sequence are within dest_alph
            if (std::all_of(sequence.cbegin(alph.alphabet_size()), sequence.cend(alph.alphabet_size()),
                    [=](auto element) { return dest_alph.contains(alph[element]); })) {
                // We have to unpack sequence to repack it later with new alphabet
                ProtoSequence<INTERNAL, STRING> unpacked =
                        internal::reserve_space_for_unpacked<INTERNAL, INTERNAL, STRING>(sequence);
                internal::unpack_common(sequence, unpacked, alph);
                return unpacked;
            } else {
                return ProtoSequence<INTERNAL, STRING>{};
            }
        }
    }

    template<InternalType INTERNAL>
    Sq<INTERNAL> remove_ambiguous(const Sq<INTERNAL> &sq,
                                  const bool by_letter) {
        const Alphabet &alph = sq.alphabet();
        SqType type;
        switch (sq.type()) {
            // TODO: should we include _BSC types as well and simply return sq?
            case AMI_EXT:
                type = AMI_BSC;
                break;
            case DNA_EXT:
                type = DNA_BSC;
                break;
            case RNA_EXT:
                type = RNA_BSC;
                break;
            default:
                throw std::invalid_argument("sq object must have extended alphabet type");
        }
        const Alphabet dest_alph = Alphabet(type, alph.NA_letter());
        Sq<INTERNAL> ret(sq.length(), dest_alph);

        for (LenSq i = 0; i < sq.length(); ++i) {
            ProtoSequence<INTERNAL, STRING> unpacked =
                    remove_ambiguous<INTERNAL>(sq[i].get(), alph, dest_alph, by_letter);
            Sequence<INTERNAL> repacked =
                    internal::reserve_space_for_packed<INTERNAL>(unpacked.length(), dest_alph.alphabet_size());
            internal::pack<INTERNAL, STRING, INTERNAL, true>(unpacked, repacked, dest_alph);
            ret[i] = repacked;
        }
        return ret;
    }
}
