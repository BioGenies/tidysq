#pragma once

#include <tidysq/Alphabet.h>
#include <tidysq/constants/ambiguous_maps.h>
#include <tidysq/Sequence.h>
#include <tidysq.h>

namespace tidysq::internal {
    template<typename INTERNAL>
    class MotifFrame;

    class Motif {
        const Alphabet &alph_;
        std::string sought_;
        std::list<std::list<unsigned short>> content_;
        bool from_start_ = false;
        bool until_end_ = false;

    private:
        [[nodiscard]] inline std::list<unsigned short> match_value(const char &letter) {
            std::list<unsigned short> ret{};
            std::list<char> meanings{};
            std::unordered_map<ElementStringSimple, std::list<ElementStringSimple>> map{};

            // Assigns mapping corresponding to sq class
            switch (alph_.type()) {
                case AMI_BSC:
                case AMI_EXT:
                    map = constants::AMBIGUOUS_AMINO_MAP;
                    break;
                case DNA_BSC:
                case DNA_EXT:
                    map = constants::AMBIGUOUS_DNA_MAP;
                    break;
                case RNA_BSC:
                case RNA_EXT:
                    map = constants::AMBIGUOUS_RNA_MAP;
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
                alph_(alph), sought_(motif) {
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

        inline static std::list<internal::Motif> convert_motifs(const std::vector<std::string>& motifs,
                                                                const Alphabet& alph) {
            std::list<internal::Motif> ret{};
            for (const auto& motif : motifs) {
                ret.emplace_back(motif, alph);
            }
            return ret;
        }

    public:
        [[nodiscard]] inline bool empty() const {
            return content_.empty();
        }

        [[nodiscard]] inline R_xlen_t size() const {
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
        [[nodiscard]] bool aligns_with(typename tidysq::Sequence<INTERNAL>::const_iterator sequence_it,
                                       const typename tidysq::Sequence<INTERNAL>::const_iterator &iterator_end) const {
            auto motif_it = begin();
            while (sequence_it <= iterator_end && std::any_of(
                    motif_it->begin(), motif_it->end(), [&](const unsigned short &possible_letter) {
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
        void locate(const Sequence<INTERNAL> &sequence,
                    typename tidysq::Sequence<INTERNAL>::const_iterator sequence_it,
                    const typename tidysq::Sequence<INTERNAL>::const_iterator &iterator_end,
                    const std::string &name,
                    MotifFrame<INTERNAL> &ret) const {
            auto motif_it = begin();
            while (sequence_it <= iterator_end && std::any_of(
                    motif_it->begin(), motif_it->end(), [&](const LetterValue & possible_letter) {
                        return *sequence_it == possible_letter;
                    })) {
                ++motif_it;
                ++sequence_it;
                // Success is only whenever we arrive at the end of the motif before the end of the sequence
                // or before motif stops corresponding to sequence
                if (motif_it == end()) {
                    // TODO: append located motif and other info to ret
                    std::vector<LenSq> indices(content_.size());
                    std::generate(indices.rbegin(), indices.rend(), [=]() mutable {
                        return (--sequence_it).index();
                    });
                    Sequence<INTERNAL> found_sequence = bite(sequence, alph_.alphabet_size(), indices);
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
            if (sequence.original_length() >= size()) {
                // Lot of ^ and $ handling mostly
                if (from_start_) {
                    if (until_end_) {
                        contains_motif = (sequence.original_length() == size()) &&
                                         aligns_with<INTERNAL>(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()));
                    } else {
                        contains_motif = aligns_with<INTERNAL>(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()));
                    }
                } else if (until_end_) {
                    contains_motif = aligns_with<INTERNAL>(sequence.cend(alph_.alphabet_size()) - size(), sequence.cend(alph_.alphabet_size()));
                } else {
                    // Basic case below (without ^ or $)
                    typename tidysq::Sequence<INTERNAL>::const_iterator it = sequence.cbegin(alph_.alphabet_size());
                    // Stop when motif no longer fits in what little part of sequence is left or we already
                    // know that there is a motif here
                    while (!contains_motif && it <= sequence.cend(alph_.alphabet_size()) - size()) {
                        contains_motif = aligns_with<INTERNAL>(it, sequence.cend(alph_.alphabet_size()));
                        ++it;
                    }
                }
            }
            return contains_motif;
        }

        template<typename INTERNAL>
        MotifFrame<INTERNAL> find_in(const Sequence<INTERNAL> &sequence,
                                     const std::string &name) const {
            auto ret = MotifFrame<INTERNAL>(alph_);
            // Don't run checks if motif is longer than sequence
            if (sequence.original_length() >= size()) {
                // Lot of ^ and $ handling mostly
                if (from_start_) {
                    if (!until_end_ || sequence.original_length() == size()) {
                        locate(sequence, sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()), name, ret);
                    }
                } else if (until_end_) {
                    locate(sequence, sequence.cend(alph_.alphabet_size()) - size(), sequence.cend(alph_.alphabet_size()), name, ret);
                } else {
                    // Basic case below (without ^ or $)
                    typename tidysq::Sequence<INTERNAL>::const_iterator it = sequence.cbegin(alph_.alphabet_size());
                    while (it <= sequence.cend(alph_.alphabet_size()) - size()) {
                        locate(sequence, it, sequence.cend(alph_.alphabet_size()), name, ret);
                        ++it;
                    }
                }
            }
            return ret;
        }
    };
}