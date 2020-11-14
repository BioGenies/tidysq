#pragma once

#include <functional>

#include "tidysq/Sq.h"
#include "tidysq/internal/pack.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> remove_on_condition(const Sequence<INTERNAL> &sequence,
                                           const Alphabet &alph,
                                           const Alphabet &dest_alph,
                                           const std::function<bool(LetterValue)> &condition,
                                           const bool by_letter) {
        typedef typename TypeMapper<INTERNAL, STRINGS_PT>::ProtoSequenceContentType ContentType;
        // TODO: extract functions out of these ifs, maybe?
        if (by_letter) {
            ContentType selected_letters;
            for (auto it = sequence.cbegin(alph.alphabet_size()); it != sequence.cend(alph.alphabet_size()); ++it) {
                if (condition(*it)) {
                    selected_letters.push_back(alph[*it]);
                }
            }
            ProtoSequence<INTERNAL, STRINGS_PT> unpacked{selected_letters};
            Sequence<INTERNAL> repacked =
                    util::reserve_space_for_packed<INTERNAL>(unpacked.length(), dest_alph.alphabet_size());
            internal::pack<INTERNAL, STRINGS_PT, INTERNAL, true>(unpacked, repacked, dest_alph);
            return repacked;
        } else {
            if (std::all_of(sequence.cbegin(alph.alphabet_size()), sequence.cend(alph.alphabet_size()),
                            [=](const LetterValue element) { return condition(element); })) {
                if (alph != dest_alph) {
                    ProtoSequence<INTERNAL, STRINGS_PT> unpacked =
                            util::reserve_space_for_unpacked<INTERNAL, INTERNAL, STRINGS_PT>(sequence);
                    internal::unpack_common(sequence, unpacked, alph);
                    Sequence<INTERNAL> repacked =
                            util::reserve_space_for_packed<INTERNAL>(unpacked.length(), dest_alph.alphabet_size());
                    internal::pack<INTERNAL, STRINGS_PT, INTERNAL, true>(unpacked, repacked, dest_alph);
                    return repacked;
                } else {
                    return sequence;
                }
            } else {
                return Sequence<INTERNAL>{};
            }
        }
    }

    template<typename INTERNAL>
    Sequence<INTERNAL> remove_NA(const Sequence<INTERNAL> &sequence,
                                 const Alphabet &alph,
                                 const bool by_letter) {
        return remove_on_condition<INTERNAL>(sequence, alph, alph, [=](LetterValue value) {
            return alph.NA_value() != value;
        }, by_letter);
    }

    template<typename INTERNAL>
    Sequence<INTERNAL> remove_ambiguous(const Sequence<INTERNAL> &sequence,
                                        const Alphabet &alph,
                                        const Alphabet &dest_alph,
                                        const bool by_letter) {
        return remove_on_condition<INTERNAL>(sequence, alph, dest_alph, [=](LetterValue value) {
            return dest_alph.contains(alph[value]) || alph.NA_value() == value;
        }, by_letter);
    }

    template<typename INTERNAL>
    Sq<INTERNAL> remove_NA(const Sq<INTERNAL> &sq,
                           const bool by_letter) {
        const Alphabet &alph = sq.alphabet();
        Sq<INTERNAL> ret(sq.length(), alph);

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = remove_NA<INTERNAL>(sq[i].get(), alph, by_letter);
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> remove_ambiguous(const Sq<INTERNAL> &sq,
                                  const bool by_letter) {
        SqType type;
        switch (sq.type()) {
            case AMI_BSC:
                return sq;
            case AMI_EXT:
                type = AMI_BSC;
                break;
            case DNA_BSC:
                return sq;
            case DNA_EXT:
                type = DNA_BSC;
                break;
            case RNA_BSC:
                return sq;
            case RNA_EXT:
                type = RNA_BSC;
                break;
            default:
                throw std::invalid_argument("sq object must have extended alphabet type");
        }
        const Alphabet &alph = sq.alphabet();
        const Alphabet dest_alph = Alphabet(type, alph.NA_letter());
        Sq<INTERNAL> ret(sq.length(), dest_alph);

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = remove_ambiguous<INTERNAL>(sq[i].get(), alph, dest_alph, by_letter);
        }
        return ret;
    }
}
