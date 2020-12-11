#pragma once 

#include <utility>

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/internal/Motif.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN>
        class OperationHas : public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                                                            std::vector<bool>, bool> {
            const std::list<internal::Motif> motif_list_;

        public:
            OperationHas(const Alphabet &alphabet,
                         const std::vector<std::string> &motifs) :
                    motif_list_(internal::Motif::convert_motifs(motifs, alphabet)) {
                // TODO: issue #61
                if (!alphabet.is_simple())
                    throw std::invalid_argument("For now, %has% is supported only for simple letter alphabets");
            };

            std::vector<bool> initialize_vector_out(
                    const Sq<INTERNAL_IN> &vector_in, const LenSq from, const LenSq to) override {
                return std::vector<bool>(to - from);
            };


            bool initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return false;
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence, bool &present) override {
                present = std::all_of(motif_list_.begin(), motif_list_.end(), [&](const internal::Motif& motif) {
                    return motif.appears_in<INTERNAL_IN>(sequence);
                });
            }

            inline bool operator() (const Sequence<INTERNAL_IN> &sequence) override {
                bool present = initialize_element_out(sequence);
                operator()(sequence, present);
                return present;
            }
        };
    }

    template<typename INTERNAL_IN>
    std::vector<bool> has(const Sq<INTERNAL_IN> &sq,
                          const std::vector<std::string> &motifs) {
        return sqapply(sq, ops::OperationHas<INTERNAL_IN>(sq.alphabet(), motifs));
    }

    template<typename INTERNAL_IN>
    bool has(const Sequence<INTERNAL_IN> &sequence,
             const Alphabet &alphabet,
             const std::vector<std::string> &motifs) {
        return ops::OperationHas<INTERNAL_IN>(alphabet, motifs)(sequence);
    }
}