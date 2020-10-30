#ifndef TIDYSQ_HAS_H
#define TIDYSQ_HAS_H

#include "tidysq/types/Sq.h"

namespace tidysq {
    namespace internal {
        class Motif;
    }

    std::list<internal::Motif> convert_motifs(const std::vector<std::string>& motifs,
                                              const Alphabet& alph);

    namespace internal {
        class Motif {
            const Alphabet &alph_;
            std::list<std::list<LetterValue>> content_;
            bool from_start_ = false;
            bool until_end_ = false;

            Motif(const std::string& motif, const Alphabet& alph) :
                    alph_(alph) {
                content_ = {};
                for (auto it = motif.begin(); it != motif.end(); ++it) {
                    if (*it == '^') {
                        if (it == motif.begin()) {
                            from_start_ = true;
                            continue;
                        }
                        else throw std::invalid_argument("'^' cannot appear anywhere other than at the beginning of motif");
                    }
                    if (*it == '$') {
                        if (it == motif.end() - 1) {
                            until_end_ = true;
                            continue;
                        }
                        else throw std::invalid_argument("'$' cannot appear anywhere other than at the end of motif");
                    }
                    content_.push_back({alph.match_value(*it)});
                }
            }

        public:
            [[nodiscard]] inline bool empty() const {
                return content_.empty();
            }

            [[nodiscard]] inline LenSq length() const {
                return content_.size();
            }

            [[nodiscard]] inline auto begin() const {
                return content_.begin();
            }

            [[nodiscard]] inline auto end() const {
                return content_.end();
            }

        private:
            // sequence_it is passed as copy, because we want a new iterator that starts from that point
            template<InternalType INTERNAL>
            [[nodiscard]] bool aligns_with(SequenceIterator<INTERNAL> sequence_it,
                                           const SequenceIterator<INTERNAL>& iterator_end) const {
                auto motif_it = begin();
                while (sequence_it <= iterator_end && std::any_of(
                        motif_it->begin(), motif_it->end(), [=](const LetterValue& possible_letter) {
                            return *sequence_it == possible_letter;
                        })) {
                    ++motif_it;
                    ++sequence_it;
                    if (motif_it == end()) {
                        return true;
                    }
                }
                return false;
            }

        public:
            template<InternalType INTERNAL>
            [[nodiscard]] bool appears_in(const Sequence<INTERNAL>& sequence) const {
                SequenceIterator<INTERNAL> it = sequence.begin(alph_);
                bool contains_motif = empty();
                if (sequence.originalLength() >= length()) {
                    if (from_start_) {
                        contains_motif = aligns_with(it, sequence.end(alph_));
                    } else {
                        while (!contains_motif && it < sequence.end(alph_) - length()) {
                            contains_motif = aligns_with(it, sequence.end(alph_));
                            ++it;
                        }
                    }
                }
                return contains_motif;
            }

            friend std::list<internal::Motif> tidysq::convert_motifs(const std::vector<std::string>& motifs,
                    const Alphabet& alph);
        };
    }

    std::list<internal::Motif> convert_motifs(const std::vector<std::string>& motifs,
                                              const Alphabet& alph) {
        std::list<internal::Motif> ret{};
        for (const auto& motif : motifs) {
            ret.push_back(internal::Motif(motif, alph));
        }
        return ret;
    }

    template<InternalType INTERNAL>
    Rcpp::LogicalVector has(const Sq<INTERNAL> &sq, const std::vector<std::string>& motifs){
        using internal::Motif;

        const Alphabet& alph = sq.alphabet();
        Rcpp::LogicalVector ret(sq.length());

        // TODO: implement possibility of reading motifs for multiletter alphabets
        if (!alph.is_simple()) throw std::exception();

        const std::list<Motif> motif_list = convert_motifs(motifs, alph);
        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = std::all_of(motif_list.begin(), motif_list.end(), [=](const Motif& motif) {
                return motif.appears_in<INTERNAL>(sq[i]);
            });
        }
        // Steps to take:
        // 1. handle ^ and $ signs
        // 2. convert motifs to a list of lists of letters
        //      (if you go for a class, place it within tidysq-internal namespace)
        //      (each element of inner lists is a set of acceptable letters at this position - stored as bits)
        //      (store one motif as linked list of arrays)
        // 3. iterate over sequences and motifs, trying to align them
        //      (probably for each vector in a list within the list use any_of(letter_vec...))
        return ret;
    }
}

#endif //TIDYSQ_HAS_H
