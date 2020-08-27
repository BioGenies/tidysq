#ifndef TIDYSQ_ALPHABET_H
#define TIDYSQ_ALPHABET_H

#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>

#include "general.h"
#include "../util/common.h"
#include "../util/transition.h"

namespace tidysq {
    class Alphabet {
    private:
        const std::vector<std::string> letters_;
        const std::string NALetter_;
        const letvalue NAValue_;
        const sizealph alphabetSize_;
        const bool simple_;

        void checkLetters() const {
            for (const auto &letter : letters_) {
                if (letter.empty()) {
                    throw std::invalid_argument(R"(All letters in the alphabet should have length greater than 0!)");
                }
            }
        }

        void checkNALetter() const {
            if (NALetter_.empty()) {
                throw std::invalid_argument(R"("na_character" should not be empty!)");
            }
        }

        [[nodiscard]] sizealph calculateAlphabetSize() const {
            return ceil(log2((double) letters_.size() + 2));
        }

        [[nodiscard]] letvalue calculateNAValue() const {
            return pow(2, alphabetSize_) - 1;
        }

        [[nodiscard]] bool calculateSimple() const {
            for (const auto &letter : letters_) {
                if (letter.size() > 1) {
                    return false;
                }
            }
            return true;
        }

    public:
        Alphabet(std::vector<std::string> letters, std::string NALetter) :
                letters_(std::move(letters)),
                NALetter_(std::move(NALetter)),
                alphabetSize_(calculateAlphabetSize()),
                NAValue_(calculateNAValue()),
                simple_(calculateSimple()) {
            checkLetters();
            checkNALetter();
        }

        explicit Alphabet(const Rcpp::StringVector &alphabet) :
                Alphabet(util::convertStringVector<0>(alphabet),
                         util::getNACharacterAsString<0>(alphabet)) {};

        explicit Alphabet(Rcpp::List::const_AttributeProxy alphabet) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(alphabet)) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) = default;

        explicit operator Rcpp::StringVector() {
            Rcpp::StringVector ret(util::convertStringVector<0>(letters_));
            ret.attr("na_character") = NALetter_;
            return ret;
        }

        [[nodiscard]] inline letvalue length() const {
            return letters_.size();
        }

        inline const std::string &operator[](letvalue index) const {
            return letters_[index];
        }

        [[nodiscard]] inline const letvalue &NAValue() const {
            return NAValue_;
        };

        [[nodiscard]] inline const std::string &NALetter() const {
            return NALetter_;
        }

        [[nodiscard]] inline const sizealph &alphabetSize() const {
            return alphabetSize_;
        }

        [[nodiscard]] inline bool isSimple() const {
            return simple_;
        }
    };
}

#endif //TIDYSQ_ALPHABET_H
