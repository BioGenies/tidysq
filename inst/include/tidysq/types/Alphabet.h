#ifndef TIDYSQ_ALPHABET_H
#define TIDYSQ_ALPHABET_H

#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>
#include <variant>

#include "tidysq/types/general.h"
#include "tidysq/types/sq_types.h"
#include "tidysq/util/common.h"

namespace tidysq {
    class Alphabet {
        const std::vector<Letter> letters_;
        const Letter NA_letter_;
        const AlphSize alphabet_size_;
        const LetterValue NA_value_;
        const SqType type_;
        const bool is_simple_;
        const std::vector<SimpleLetter> simple_letters_;
        const SimpleLetter simple_NA_letter_;

        void check_letters() const {
            for (auto &letter : letters_) {
                if (letters_.empty())
                    throw std::invalid_argument("each \"letter\" has to have at least one character!");
            }
        }

        void check_NA_letter() const {
            if (NA_letter_.empty())
                throw std::invalid_argument("\"NA_letter\" has to have at least one character!");
        }

        [[nodiscard]] AlphSize calculate_alphabet_size() const {
            return ceil(log2((double) letters_.size() + 1));
        }

        [[nodiscard]] LetterValue calculate_NA_value() const {
            return pow(2, alphabet_size_) - 1;
        }

        [[nodiscard]] bool calculate_is_simple() const {
            return NA_letter_.size() == 1 &&
                std::all_of(letters_.begin(), letters_.end(), [](const Letter& letter){ return letter.size() == 1; });
        }

        [[nodiscard]] std::vector<char> create_simple_letters() const {
            if (!is_simple_) return {};
            auto ret = std::vector<char>(letters_.size());
            for(LetterValue i = 0; i < letters_.size(); i++) {
                ret[i] = letters_[i][0];
            }
            return ret;
        }

        [[nodiscard]] char create_simple_NA_letter() const {
            return NA_letter_[0];
        }

    public:
        Alphabet(const std::vector<Letter> &letters,
                 const SqType &type,
                 const Letter &NA_letter = util::default_NA_letter()) :
                letters_(letters),
                NA_letter_(NA_letter),
                alphabet_size_(calculate_alphabet_size()),
                NA_value_(calculate_NA_value()),
                is_simple_(calculate_is_simple()),
                simple_letters_(create_simple_letters()),
                simple_NA_letter_(create_simple_NA_letter()),
                type_(type) {
            check_letters();
            check_NA_letter();
        }

        explicit Alphabet(const SqType &type,
                          const Letter &NA_letter = util::default_NA_letter()) :
                      Alphabet(util::standard_letters_for_sq_type(type),
                              type,
                              NA_letter) {};

        explicit Alphabet(const std::vector<Letter> &letters,
                          const Letter &NA_letter = util::default_NA_letter()) :
                Alphabet(util::guess_sq_type_from_letters(letters),
                         NA_letter) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) noexcept = default;

        [[nodiscard]] inline LetterValue length() const {
            return letters_.size();
        }

        inline const Letter &operator[](LetterValue index) const {
            return index == NA_value_ ? NA_letter_ : letters_[index];
        }

        inline const SimpleLetter &get_simple_letter(LetterValue index) const {
            return index == NA_value_ ? simple_NA_letter_ : simple_letters_[index];
        }

        [[nodiscard]] inline const std::vector<Letter> &letters() const {
            return letters_;
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

        inline bool operator==(const Alphabet &other) const {
            return letters_ == other.letters_;
        }

        inline bool operator!=(const Alphabet &other) const {
            return !operator==(other);
        }

        [[nodiscard]] inline LetterValue match_value(const ElementRaws &letter) const {
            if (letter < letters_.size()) return letter;
            return NA_value_;
        }

        [[nodiscard]] inline LetterValue match_value(const ElementInts &letter) const {
            if (letter < letters_.size()) return letter;
            return NA_value_;
        }

        [[nodiscard]] inline LetterValue match_value(const ElementStringSimple &letter) const {
            for(LetterValue i = 0; i < letters_.size(); i++) {
                if (letter == simple_letters_[i]) return i;
            }
            return NA_value_;
        }

        [[nodiscard]] inline LetterValue match_value(const Letter &letter) const {
            for(LetterValue i = 0; i < letters_.size(); i++) {
                if (letter == letters_[i]) return i;
            }
            return NA_value_;
        }

        friend Rcpp::StringVector export_to_R(const Alphabet &alphabet);
    };
}

#endif //TIDYSQ_ALPHABET_H
