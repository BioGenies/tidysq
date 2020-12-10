#pragma once 

#include "tidysq/internal/MotifFrame.h"
#include "tidysq/internal/Motif.h"

namespace tidysq {
    template<typename INTERNAL>
    internal::MotifFrame<INTERNAL> find_motifs(const Sq<INTERNAL> &sq,
                                                const std::vector<std::string>& names,
                                                const std::vector<std::string>& motifs) {
        const Alphabet& alph = sq.alphabet();
        // TODO: implement possibility of reading motifs for multiletter alphabets
        if (!alph.is_simple())
            throw std::invalid_argument("For now, 'find_motifs' is supported only for simple letter alphabets");

        auto motif_list = internal::Motif::convert_motifs(motifs, alph);
        internal::MotifFrame<INTERNAL> ret(alph);

        for (const internal::Motif &motif : motif_list) {
            for (LenSq i = 0; i < sq.size(); ++i) {
                ret.merge_with(motif.find_in<INTERNAL>(sq[i], names[i]));
            }
        }
        // Things to return:
        // name: name of the sequence, in which the motif was found
        // found: part of sequence that aligns with the motif
        // sought: aligned motif
        // start: beginning of alignment
        // end: end of alignment
        return ret;
    }
}