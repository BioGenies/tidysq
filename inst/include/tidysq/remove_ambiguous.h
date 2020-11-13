#pragma once

#include <functional>

#include "tidysq/exports.h"

namespace tidysq {
    template<InternalType INTERNAL>
    ProtoSequence<INTERNAL, STRINGS> remove_on_condition(const Sequence<INTERNAL> &sequence,
                                                        const Alphabet &alph,
                                                        const std::function<bool(LetterValue)> &condition,
                                                        const bool by_letter) {
        typedef typename TypeMapper<INTERNAL, STRINGS>::ProtoSequenceContentType ContentType;
        // TODO: extract functions out of these ifs, maybe?
        if (by_letter) {
            ContentType selected_letters;
            for (auto it = sequence.cbegin(alph.alphabet_size()); it != sequence.cend(alph.alphabet_size()); ++it) {
                if (condition(*it)) {
                    selected_letters.push_back(alph[*it]);
                }
            }
            return ProtoSequence<INTERNAL, STRINGS>{selected_letters};
        } else {
            if (std::all_of(sequence.cbegin(alph.alphabet_size()), sequence.cend(alph.alphabet_size()),
                            [=](const LetterValue element) { return condition(element); })) {
                // TODO: now we have to unpack sequence to repack it later with new alphabet;
                //  would be nice not to
                ProtoSequence<INTERNAL, STRINGS> unpacked =
                        internal::reserve_space_for_unpacked<INTERNAL, INTERNAL, STRINGS>(sequence);
                internal::unpack_common(sequence, unpacked, alph);
                return unpacked;
            } else {
                return ProtoSequence<INTERNAL, STRINGS>{};
            }
        }
    }

    template<InternalType INTERNAL>
    ProtoSequence<INTERNAL, STRINGS> remove_ambiguous(const Sequence<INTERNAL> &sequence,
                                                     const Alphabet &alph,
                                                     const Alphabet &dest_alph,
                                                     const bool by_letter) {
        return remove_on_condition<INTERNAL>(sequence, alph, [=](LetterValue value) {
            return dest_alph.contains(alph[value]) || alph.NA_value() == value;
        }, by_letter);
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
            ProtoSequence<INTERNAL, STRINGS> unpacked =
                    remove_ambiguous<INTERNAL>(sq[i].get(), alph, dest_alph, by_letter);
            Sequence<INTERNAL> repacked =
                    internal::reserve_space_for_packed<INTERNAL>(unpacked.length(), dest_alph.alphabet_size());
            internal::pack<INTERNAL, STRINGS, INTERNAL, true>(unpacked, repacked, dest_alph);
            ret[i] = repacked;
        }
        return ret;
    }
}
