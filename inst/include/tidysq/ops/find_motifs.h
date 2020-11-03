#ifndef TIDYSQ_FIND_MOTIFS_H
#define TIDYSQ_FIND_MOTIFS_H

#include "tidysq/types/Sq.h"
#include <map>

namespace tidysq {
    namespace internal {
        class Motif;
    }

    typedef std::unordered_map<ElementStringSimple, std::list<ElementStringSimple>> AmbiguousDict;

    std::list<internal::Motif> convert_motifs(const std::vector<std::string>& motifs,
                                              const Alphabet& alph);

    AmbiguousDict ambiguousAminoMap = {
            {'B', {'B', 'D', 'N'}},
            {'J', {'J', 'I', 'L'}},
            {'Z', {'Z', 'E', 'Q'}},
            {'X', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
                          'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}}
    };

    AmbiguousDict ambiguousDNAMap = {
            {'W', {'W', 'A', 'T'}},
            {'S', {'S', 'C', 'G'}},
            {'M', {'M', 'A', 'C'}},
            {'K', {'K', 'G', 'T'}},
            {'R', {'R', 'A', 'G'}},
            {'Y', {'Y', 'C', 'T'}},
            {'B', {'B', 'S', 'K', 'Y', 'C', 'G', 'T'}},
            {'D', {'D', 'W', 'K', 'R', 'A', 'G', 'T'}},
            {'H', {'H', 'W', 'M', 'Y', 'A', 'C', 'T'}},
            {'V', {'V', 'S', 'M', 'R', 'A', 'C', 'G'}},
            {'N', {'A', 'C', 'G', 'T', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'}}
    };

    AmbiguousDict ambiguousRNAMap = {
            {'W', {'W', 'A', 'U'}},
            {'S', {'S', 'C', 'G'}},
            {'M', {'M', 'A', 'C'}},
            {'K', {'K', 'G', 'U'}},
            {'R', {'R', 'A', 'G'}},
            {'Y', {'Y', 'C', 'U'}},
            {'B', {'B', 'S', 'K', 'Y', 'C', 'G', 'U'}},
            {'D', {'D', 'W', 'K', 'R', 'A', 'G', 'U'}},
            {'H', {'H', 'W', 'M', 'Y', 'A', 'C', 'U'}},
            {'V', {'V', 'S', 'M', 'R', 'A', 'C', 'G'}},
            {'N', {'A', 'C', 'G', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'}}
    };

    namespace internal {
        class Motif {
            const Alphabet &alph_;
            std::list<std::list<LetterValue>> content_;
            bool from_start_ = false;
            bool until_end_ = false;

        private:
            [[nodiscard]] inline std::list<LetterValue> match_value(const ElementStringSimple &letter) {
                std::list<LetterValue> ret{};
                std::list<ElementStringSimple> meanings{};
                AmbiguousDict map{};

                // Assigns mapping corresponding to sq class
                switch (alph_.type()) {
                    case AMI_BSC:
                    case AMI_EXT:
                        map = ambiguousAminoMap;
                        break;
                    case DNA_BSC:
                    case DNA_EXT:
                        map = ambiguousDNAMap;
                        break;
                    case RNA_BSC:
                    case RNA_EXT:
                        map = ambiguousRNAMap;
                        break;
                    default:
                        break;
                }

                // TODO: replace with .contains once C++20 becomes widely acceptable
                // Replaces each input letter with a list of letters that are encompassed by the meaning of the input letter
                if (!map.empty() && map.count(letter) == 1) {
                    meanings = ambiguousAminoMap[letter];
                } else {
                    meanings.push_back(letter);
                }

                // Translates each meaning (a char) to LetterValue (bit representation)
                for (auto meaning : meanings) {
                    ret.push_back(alph_.match_value(meaning));
                }
                return ret;
            }

        public:
            Motif(const std::string& motif, const Alphabet& alph) :
                    alph_(alph) {
                content_ = {};
                for (auto it = motif.begin(); it != motif.end(); ++it) {
                    // In general, special handling of ^ and $ -- the only regex options implemented
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
                    // match_value returns a list of bit-packed meanings
                    content_.push_back(match_value(*it));
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
                    // Success is only whenever we arrive at the end of the motif before the end of the sequence
                    // or before motif stops corresponding to sequence
                    if (motif_it == end()) {
                        return true;
                    }
                }
                return false;
            }

        public:
            template<InternalType INTERNAL>
            [[nodiscard]] bool appears_in(const Sequence<INTERNAL>& sequence) const {
                bool contains_motif = empty();
                // Don't run checks if motif is longer than sequence
                if (sequence.originalLength() >= length()) {
                    // Lot of ^ and $ handling mostly
                    if (from_start_) {
                        if (until_end_) {
                            contains_motif = (sequence.originalLength() == length()) &&
                                    aligns_with(sequence.begin(alph_), sequence.end(alph_));
                        } else {
                            contains_motif = aligns_with(sequence.begin(alph_), sequence.end(alph_));
                        }
                    } else if (until_end_) {
                        contains_motif = aligns_with(sequence.end(alph_) - length(), sequence.end(alph_));
                    } else {
                        // Basic case below (without ^ or $)
                        SequenceIterator<INTERNAL> it = sequence.begin(alph_);
                        // Stop when motif no longer fits in what little part of sequence is left or we already
                        // know that there is a motif here
                        while (!contains_motif && it <= sequence.end(alph_) - length()) {
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
            ret.emplace_back(motif, alph);
        }
        return ret;
    }

    template<InternalType INTERNAL>
    Rcpp::List find_motifs(const Sq<INTERNAL> &sq,
                           const std::vector<std::string>& names,
                           const std::vector<std::string>& motifs) {
        using internal::Motif;

        const Alphabet& alph = sq.alphabet();
        // TODO: implement possibility of reading motifs for multiletter alphabets
        if (!alph.is_simple())
            throw std::invalid_argument("For now, %has% is supported only for simple letter alphabets");

        const std::list<Motif> motif_list = convert_motifs(motifs, alph);
        
        return Rcpp::List::create();
    }

    template<InternalType INTERNAL>
    Rcpp::LogicalVector has(const Sq<INTERNAL> &sq, const std::vector<std::string>& motifs) {
        using internal::Motif;

        const Alphabet& alph = sq.alphabet();
        // TODO: implement possibility of reading motifs for multiletter alphabets
        if (!alph.is_simple())
            throw std::invalid_argument("For now, %has% is supported only for simple letter alphabets");

        const std::list<Motif> motif_list = convert_motifs(motifs, alph);
        Rcpp::LogicalVector ret(sq.length());

        for (LenSq i = 0; i < sq.length(); ++i) {
            // all_of guarantees early stopping if any motif is not present
            ret[i] = std::all_of(motif_list.begin(), motif_list.end(), [=](const Motif& motif) {
                return motif.appears_in<INTERNAL>(sq[i]);
            });
        }
        // Steps to take:
        // 1. add support for multicharacter alphabet
        return ret;
    }
}

#endif //TIDYSQ_FIND_MOTIFS_H
