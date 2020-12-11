#pragma once

#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/sq-types.h"
#include "tidysq/util/transform-common.h"
#include "tidysq/constants/standard_letters.h"

namespace tidysq {
    template<typename INTERNAL>
    class Sq;

    class Alphabet {
        bool ignore_case_;
        std::unordered_map<LetterValue, const Letter> value_to_letter_;
        Letter NA_letter_;
        AlphSize alphabet_size_;
        LetterValue NA_value_;
        bool is_simple_;
        std::unordered_map<LetterValue, SimpleLetter> value_to_simple_letter_;
        SimpleLetter NA_simple_letter_;
        std::unordered_map<Letter, LetterValue> letter_to_value_;
        std::unordered_map<SimpleLetter, LetterValue> simple_letter_to_value_;
        SqType type_;

        [[nodiscard]] inline AlphSize calculate_alphabet_size() {
            return static_cast<AlphSize>(ceil(log2((double) value_to_letter_.size() + 1)));
        }

        [[nodiscard]] inline LetterValue calculate_NA_value() const {
            return static_cast<LetterValue>(pow(2, alphabet_size_) - 1);
        }

        [[nodiscard]] inline bool calculate_is_simple() const {
            return NA_letter_.size() == 1 &&
                std::all_of(value_to_letter_.begin(), value_to_letter_.end(), [](const auto& pair){ return pair.second.size() == 1; });
        }

        [[nodiscard]] SimpleLetter prepare_simple_NA_letter() const {
            return NA_letter_[0];
        }

        static std::unordered_map<LetterValue, const Letter> prepare_value_to_letter(const std::vector<Letter> &letters) {
            std::unordered_map<LetterValue, const Letter> ret{};
            for (int i = 0; i < letters.size(); i++) {
                if (letters[i].empty())
                    throw std::invalid_argument("each \"letter\" has to have at least one character!");
                ret.insert({i, letters[i]});
            }
            return ret;
        }

        std::unordered_map<Letter, LetterValue> prepare_letter_to_value() {
            std::unordered_map<Letter, LetterValue> ret{};
            if (ignore_case_) {
                if (!is_simple_)
                    throw std::invalid_argument("\"ignore_case\" cannot be used with non-simple alphabet");
                for (const auto &pair : value_to_letter_) {
                    ret.insert({pair.second, pair.first});
                    if (tolower(pair.second[0]) != pair.second[0])
                        ret.insert({std::string{(char) tolower(pair.second[0])}, pair.first});
                }
            } else {
                for (const auto &pair : value_to_letter_) {
                    ret.insert({pair.second, pair.first});
                }
            }
            return ret;
        }

        [[nodiscard]] std::unordered_map<SimpleLetter, LetterValue> prepare_simple_letter_to_value() const {
            if (!is_simple_) return {};
            std::unordered_map<SimpleLetter, LetterValue> ret{};
            for(const auto &pair : letter_to_value_) {
                ret.insert({pair.first[0], pair.second});
            }
            return ret;
        }

        static inline Letter prepare_NA_letter(const Letter &NA_letter) {
            if (NA_letter.empty())
                throw std::invalid_argument("\"NA_letter\" has to have at least one character!");
            return NA_letter;
        }

        [[nodiscard]] std::unordered_map<LetterValue, SimpleLetter> prepare_value_to_simple_letter() const {
            if (!is_simple_) return {};
            std::unordered_map<LetterValue, SimpleLetter> ret{};
            for (const auto &pair : value_to_letter_) {
                ret.insert({pair.first, pair.second[0]});
            }
            return ret;
        }

    public:
        typedef typename std::unordered_map<LetterValue, const Letter>::const_iterator const_iterator;

        Alphabet(const std::vector<Letter> &letters,
                 const SqType &type,
                 const Letter &NA_letter = constants::DEFAULT_NA_LETTER,
                 const bool ignore_case = constants::DEFAULT_IGNORE_CASE) :
                ignore_case_(ignore_case),
                value_to_letter_(prepare_value_to_letter(letters)),
                NA_letter_(prepare_NA_letter(NA_letter)),
                alphabet_size_(calculate_alphabet_size()),
                NA_value_(calculate_NA_value()),
                is_simple_(calculate_is_simple()),
                value_to_simple_letter_(prepare_value_to_simple_letter()),
                NA_simple_letter_(prepare_simple_NA_letter()),
                letter_to_value_(prepare_letter_to_value()),
                simple_letter_to_value_(prepare_simple_letter_to_value()),
                type_(type) {}

        explicit Alphabet(const SqType &type,
                          const Letter &NA_letter = constants::DEFAULT_NA_LETTER,
                          const bool ignore_case = constants::DEFAULT_IGNORE_CASE) :
                      Alphabet(util::standard_letters_for_sq_type(type),
                              type,
                              NA_letter,
                              ignore_case) {};

        explicit Alphabet(const std::vector<Letter> &letters,
                          const Letter &NA_letter = constants::DEFAULT_NA_LETTER,
                          const bool ignore_case = constants::DEFAULT_IGNORE_CASE) :
                Alphabet(util::has_standard_alphabet(util::guess_sq_type_from_letters(letters)) ?
                util::standard_letters_for_sq_type(util::guess_sq_type_from_letters(letters)) : letters,
                         util::guess_sq_type_from_letters(letters),
                         NA_letter,
                         ignore_case) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) noexcept = default;

        Alphabet& operator=(const Alphabet &other) = default;

        [[nodiscard]] inline LetterValue size() const {
            return static_cast<LetterValue>(value_to_letter_.size());
        }

        inline const Letter &operator[](LetterValue index) const {
            return index == NA_value_ ? NA_letter_ : value_to_letter_.at(index);
        }

        inline const SimpleLetter &get_simple_letter(LetterValue index) const {
            return index == NA_value_ ? NA_simple_letter_ : value_to_simple_letter_.at(index);
        }

        [[nodiscard]] inline const LetterValue &NA_value() const {
            return NA_value_;
        }

        [[nodiscard]] inline const Letter &NA_letter() const {
            return NA_letter_;
        }

        [[nodiscard]] inline const SqType &type() const {
            return type_;
        }

        [[nodiscard]] inline const AlphSize &alphabet_size() const {
            return alphabet_size_;
        }

        [[nodiscard]] inline bool is_simple() const {
            return is_simple_;
        }

        [[nodiscard]] inline bool ignores_case() const {
            return ignore_case_;
        }

        [[nodiscard]] inline const_iterator begin() const {
            return value_to_letter_.begin();
        }

        [[nodiscard]] inline const_iterator end() const {
            return value_to_letter_.end();
        }

        [[nodiscard]] inline const_iterator cbegin() const {
            return value_to_letter_.cbegin();
        }

        [[nodiscard]] inline const_iterator cend() const {
            return value_to_letter_.cend();
        }

        inline bool operator==(const Alphabet &other) const {
            return value_to_letter_ == other.value_to_letter_ &&
                     NA_letter_ == other.NA_letter_;
        }

        inline bool operator!=(const Alphabet &other) const {
            return !operator==(other);
        }

        [[nodiscard]] inline LetterValue match_value(const ElementRaws &letter) const {
            if (letter < value_to_letter_.size()) return letter;
            return NA_value_;
        }

        [[maybe_unused]] [[nodiscard]] inline LetterValue match_value(const ElementInts &letter) const {
            if (letter < value_to_letter_.size()) return letter;
            return NA_value_;
        }

        [[nodiscard]] inline LetterValue match_value(const ElementStringSimple &letter) const {
            try {
                return simple_letter_to_value_.at(letter);
            } catch (const std::out_of_range &e) {
                return NA_value_;
            }
        }

        [[nodiscard]] inline LetterValue match_value(const Letter &letter) const {
            try {
                return letter_to_value_.at(letter);
            } catch (const std::out_of_range &e) {
                return NA_value_;
            }
        }

        [[nodiscard]] inline bool contains(const Letter &letter) const {
            return std::any_of(cbegin(), cend(), [=](const std::pair<LetterValue, const Letter> &entry) {
                return letter == entry.second;
            });
        }

        friend Rcpp::StringVector export_to_R(const Alphabet &alphabet);
    };
}
