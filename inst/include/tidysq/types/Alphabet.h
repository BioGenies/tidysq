#ifndef TIDYSQ_ALPHABET_H
#define TIDYSQ_ALPHABET_H

#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>

#include "tidysq/types/general.h"
#include "tidysq/util/common.h"
#include "tidysq/util/transition.h"

namespace tidysq {
    class Alphabet {
    private:
        const std::vector<std::string> letters_;
        const std::string NALetter_;
        const LetValue NAValue_;
        const AlphSize alphabetSize_;
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
                throw std::invalid_argument(R"("na_letter" should not be empty!)");
            }
        }

        static std::vector<std::string> convertSetToVector(const std::set<char> &letters) {
            std::vector<std::string> ret(letters.size());
            auto iterator = letters.begin();
            while (iterator != letters.end()) {
                ret.emplace_back(std::string({*iterator}));
                iterator++;
            }
            return ret;
        }

        [[nodiscard]] AlphSize calculateAlphabetSize() const {
            return ceil(log2((double) letters_.size() + 2));
        }

        [[nodiscard]] LetValue calculateNAValue() const {
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
                Alphabet(util::convertStringVector(alphabet),
                         util::getNACharacterAsString(alphabet)) {};

        explicit Alphabet(Rcpp::List::const_AttributeProxy alphabet) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(alphabet)) {};

        Alphabet(const std::set<char> &letters, std::string NALetter) :
                Alphabet(convertSetToVector(letters), std::move(NALetter)) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) = default;

        explicit operator Rcpp::StringVector() {
            Rcpp::StringVector ret(util::convertStringVector(letters_));
            ret.attr("na_letter") = NALetter_;
            return ret;
        }

        [[nodiscard]] inline LetValue length() const {
            return letters_.size();
        }

        inline const std::string &operator[](LetValue index) const {
            return letters_[index];
        }

        [[nodiscard]] inline const LetValue &NAValue() const {
            return NAValue_;
        }

        [[nodiscard]] inline const std::string &NALetter() const {
            return NALetter_;
        }

        [[nodiscard]] inline const AlphSize &alphabetSize() const {
            return alphabetSize_;
        }

        [[nodiscard]] inline bool isSimple() const {
            return simple_;
        }

        bool operator==(const Alphabet &another) {
            if (letters_ == another.letters_) return true;
            return false;
        }
    };
}

#endif //TIDYSQ_ALPHABET_H
