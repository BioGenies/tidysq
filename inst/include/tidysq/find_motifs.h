#pragma once
#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err58-cpp"

#include "tidysq/Sq.h"
#include "tidysq/bite.h"
#include <map>
#include <algorithm>

namespace tidysq {
    namespace internal {
        template<typename INTERNAL>
        class FoundMotifs;
        class Motif;
    }

    typedef std::unordered_map<ElementStringSimple, std::list<ElementStringSimple>> AmbiguousDict;

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
        template<typename INTERNAL>
        class FoundMotifs {
            std::list<std::string> names_{};
            Sq<INTERNAL> found_;
            std::list<std::string> sought_{};
            std::list<LenSq> start_{};
            std::list<LenSq> end_{};

        public:
            explicit FoundMotifs(const Sq<INTERNAL> &sq) :
                    found_(Sq<INTERNAL>(sq.alphabet())) {}

            inline void append(const std::string &name,
                               const Sequence<INTERNAL> &found,
                               const std::string &sought,
                               const LenSq start,
                               const LenSq end) {
                names_.push_back(name);
                found_.push_back(found);
                sought_.push_back(sought);
                start_.push_back(start);
                end_.push_back(end);
            }

            friend Rcpp::List export_to_R(const internal::FoundMotifs<RCPP_IT> &found_motifs);
        };

        Rcpp::List export_to_R(const internal::FoundMotifs<RCPP_IT> &found_motifs) {

            return Rcpp::List::create(
                    Rcpp::Named("names", found_motifs.names_),
                    Rcpp::Named("found", export_to_R(found_motifs.found_)),
                    Rcpp::Named("sought", found_motifs.sought_),
                    Rcpp::Named("start", Rcpp::IntegerVector(Rcpp::wrap(found_motifs.start_)) + 1),
                    Rcpp::Named("end", Rcpp::IntegerVector(Rcpp::wrap(found_motifs.end_)) + 1));
        }

        class Motif {
            const Alphabet &alph_;
            std::string sought_;
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
                    meanings = map[letter];
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
            Motif(const std::string &motif, const Alphabet &alph) :
                    sought_(motif), alph_(alph) {
                content_ = {};
                for (auto it = motif.begin(); it != motif.end(); ++it) {
                    // In general, special handling of ^ and $ -- the only regex options implemented
                    if (*it == '^') {
                        if (it == motif.begin()) {
                            from_start_ = true;
                            continue;
                        } else
                            throw std::invalid_argument(
                                    "'^' cannot appear anywhere other than at the beginning of motif");
                    }
                    if (*it == '$') {
                        if (it == motif.end() - 1) {
                            until_end_ = true;
                            continue;
                        } else throw std::invalid_argument("'$' cannot appear anywhere other than at the end of motif");
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
            template<typename INTERNAL>
            [[nodiscard]] bool aligns_with(typename Sequence<INTERNAL>::const_iterator sequence_it,
                                           const typename Sequence<INTERNAL>::const_iterator &iterator_end) const {
                auto motif_it = begin();
                while (sequence_it <= iterator_end && std::any_of(
                        motif_it->begin(), motif_it->end(), [=](const LetterValue &possible_letter) {
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

            // sequence_it is passed as copy, because we want a new iterator that starts from that point
            template<typename INTERNAL>
            void locate(typename Sequence<INTERNAL>::const_iterator sequence_it,
                        const typename Sequence<INTERNAL>::const_iterator &iterator_end,
                        const std::string &name,
                        internal::FoundMotifs<INTERNAL> &ret) const {
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
                        // TODO: append located motif and other info to ret
                        std::vector<int> indices(content_.size());
                        std::generate(indices.rbegin(), indices.rend(), [=]() mutable {
                            return (--sequence_it).index();
                        });
                        Sequence<INTERNAL> found_sequence = bite<INTERNAL>(sequence_it, indices);
                        ret.append(name, found_sequence, sought_, sequence_it.index() - content_.size(), sequence_it.index() - 1);
                        return;
                    }
                }
            }

        public:
            template<typename INTERNAL>
            [[nodiscard]] bool appears_in(const Sequence<INTERNAL>& sequence) const {
                bool contains_motif = empty();
                // Don't run checks if motif is longer than sequence
                if (sequence.original_length() >= length()) {
                    // Lot of ^ and $ handling mostly
                    if (from_start_) {
                        if (until_end_) {
                            contains_motif = (sequence.original_length() == length()) &&
                                             aligns_with<INTERNAL>(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()));
                        } else {
                            contains_motif = aligns_with<INTERNAL>(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()));
                        }
                    } else if (until_end_) {
                        contains_motif = aligns_with<INTERNAL>(sequence.cend(alph_.alphabet_size()) - length(), sequence.cend(alph_.alphabet_size()));
                    } else {
                        // Basic case below (without ^ or $)
                        typename Sequence<INTERNAL>::const_iterator it = sequence.cbegin(alph_.alphabet_size());
                        // Stop when motif no longer fits in what little part of sequence is left or we already
                        // know that there is a motif here
                        while (!contains_motif && it <= sequence.cend(alph_.alphabet_size()) - length()) {
                            contains_motif = aligns_with<INTERNAL>(it, sequence.cend(alph_.alphabet_size()));
                            ++it;
                        }
                    }
                }
                return contains_motif;
            }

            template<typename INTERNAL>
            void find_in(const Sequence<INTERNAL> &sequence,
                         const std::string &name,
                         internal::FoundMotifs<INTERNAL> &ret) const {
                // Don't run checks if motif is longer than sequence
                if (sequence.original_length() >= length()) {
                    // Lot of ^ and $ handling mostly
                    if (from_start_) {
                        if (!until_end_ || sequence.original_length() == length()) {
                            locate(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()), name, ret);
                        }
                    } else if (until_end_) {
                        locate(sequence.cend(alph_.alphabet_size()) - length(), sequence.cend(alph_.alphabet_size()), name, ret);
                    } else {
                        // Basic case below (without ^ or $)
                        typename Sequence<INTERNAL>::const_iterator it = sequence.cbegin(alph_.alphabet_size());
                        while (it <= sequence.cend(alph_.alphabet_size()) - length()) {
                            locate(it, sequence.cend(alph_.alphabet_size()), name, ret);
                            ++it;
                        }
                    }
                }
            }
        };
    }

    inline std::list<internal::Motif> convert_motifs(const std::vector<std::string>& motifs,
                                                     const Alphabet& alph) {
        std::list<internal::Motif> ret{};
        for (const auto& motif : motifs) {
            ret.emplace_back(motif, alph);
        }
        return ret;
    }

    template<typename INTERNAL>
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
        return ret;
    }

    template<typename INTERNAL>
    internal::FoundMotifs<INTERNAL> find_motifs(const Sq<INTERNAL> &sq,
                                                const std::vector<std::string>& names,
                                                const std::vector<std::string>& motifs) {
        using internal::Motif;

        const Alphabet& alph = sq.alphabet();
        // TODO: implement possibility of reading motifs for multiletter alphabets
        if (!alph.is_simple())
            throw std::invalid_argument("For now, 'find_motifs' is supported only for simple letter alphabets");

        const std::list<Motif> motif_list = convert_motifs(motifs, alph);
        internal::FoundMotifs<INTERNAL> ret(sq);

        for (const Motif &motif : motif_list) {
            for (LenSq i = 0; i < sq.length(); ++i) {
                motif.find_in<INTERNAL>(sq[i], names[i], ret);
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

#pragma clang diagnostic pop