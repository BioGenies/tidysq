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

namespace tidysq {
    namespace util {
        SqType guessSqType(const std::vector<std::string> &letters);
        std::vector<std::string> getStandardLetters(const SqType &type);
        std::vector<std::string> getClassStringVector(const SqType &type);
        SqType getSqType(const Rcpp::StringVector &classVector);
        std::string getDefaultNA_letter();
        std::string getSqTypeAbbr(const SqType &type);
        SqType get_sq_type_from_abbr(const Rcpp::StringVector &type_vector);
    }

    class Alphabet {
    private:
        const std::vector<std::string> letters_;
        const std::string NA_letter_;
        const AlphSize alphabetSize_;
        const LetValue NAValue_;
        const bool simple_;
        const SqType type_;

        void checkLetters() const {
            for (const auto &letter : letters_) {
                if (letter.empty()) {
                    throw std::invalid_argument(R"(All letters in the alphabet should have length greater than 0!)");
                }
            }
        }

        void checkNA_letter() const {
            if (NA_letter_.empty()) {
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
            return ceil(log2((double) letters_.size() + 1));
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
        Alphabet(const std::vector<std::string> &letters,
                 const SqType &type,
                 const std::string &NA_letter = util::getDefaultNA_letter()) :
                letters_(letters),
                NA_letter_(NA_letter),
                alphabetSize_(calculateAlphabetSize()),
                NAValue_(calculateNAValue()),
                simple_(calculateSimple()),
                type_(type) {
            checkLetters();
            checkNA_letter();
        }

        Alphabet(const Rcpp::StringVector &letters,
                 const SqType &type,
                 const Rcpp::StringVector &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(util::convertStringVector(letters),
                         type,
                         util::getScalarStringValue(NA_letter)) {};

        Alphabet(const Rcpp::List::const_AttributeProxy &letters,
                 const SqType &type,
                 const Rcpp::StringVector &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(letters),
                         type,
                         util::getScalarStringValue(NA_letter)) {};

        Alphabet(const std::set<char> &letters,
                 const SqType &type,
                 const std::string &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(convertSetToVector(letters),
                         type,
                         NA_letter) {};

        explicit Alphabet(const std::vector<std::string> &letters,
                 const std::string &NA_letter = util::getDefaultNA_letter()) :
                 Alphabet(letters,
                          util::guessSqType(letters),
                          NA_letter) {};

        explicit Alphabet(const Rcpp::StringVector &alphabet,
                 const Rcpp::StringVector &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(util::convertStringVector(alphabet),
                         util::get_sq_type_from_abbr(alphabet.attr("type")),
                         util::getScalarStringValue(NA_letter)) {};

        explicit Alphabet(const Rcpp::List::const_AttributeProxy &alphabet,
                 const Rcpp::StringVector &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(alphabet),
                         util::get_sq_type_from_abbr(Rcpp::as<Rcpp::StringVector>(alphabet).attr("type")),
                         util::getScalarStringValue(NA_letter)) {};

        explicit Alphabet(const std::set<char> &letters,
                 const std::string &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(convertSetToVector(letters),
                         NA_letter) {};

        explicit Alphabet(const SqType &type,
                          const std::string &NA_letter = util::getDefaultNA_letter()) :
                Alphabet(util::getStandardLetters(type), type, NA_letter) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) = default;

        explicit operator Rcpp::StringVector() {
            Rcpp::StringVector ret(util::convertStringVector(letters_));
            ret.attr("type") = util::getSqTypeAbbr(type_);
            ret.attr("class") = Rcpp::StringVector{"sq_alphabet", "character", "vctrs_vctr"};
            return ret;
        }

        [[nodiscard]] inline LetValue length() const {
            return letters_.size();
        }

        inline const std::string &operator[](LetValue index) const {
            return index == NAValue_ ? NA_letter_ : letters_[index];
        }

        [[nodiscard]] inline const LetValue &NAValue() const {
            return NAValue_;
        }

        [[nodiscard]] inline const std::string &NA_letter() const {
            return NA_letter_;
        }

        [[nodiscard]] inline const SqType &type() const {
            return type_;
        }

        [[nodiscard]] inline const AlphSize &alphabetSize() const {
            return alphabetSize_;
        }

        [[nodiscard]] inline bool isSimple() const {
            return simple_;
        }

        inline bool operator==(const Alphabet &other) const {
            return letters_ == other.letters_;
        }

        inline bool operator!=(const Alphabet &other) const {
            return !operator==(other);
        }
    };

    namespace util {
        inline std::string getDefaultNA_letter() {
            return "!";
        }

        inline std::vector<std::string> getStandardLetters(const SqType &type) {
            std::vector<std::string> letters;
            switch (type) {
                case AMI:
                    letters = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                               "S", "T", "U", "V", "W", "X", "Y", "Z", "-", "*"};
                    break;
                case AMI_CLN:
                    letters = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                               "W", "Y", "-", "*"};
                    break;
                case DNA:
                    letters = {"A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                    break;
                case DNA_CLN:
                    letters = {"A", "C", "G", "T", "-"};
                    break;
                case RNA:
                    letters = {"A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                    break;
                case RNA_CLN:
                    letters = {"A", "C", "G", "U", "-"};
                    break;
                default:
                    throw std::invalid_argument("Provided type does not have a predefined standard alphabet!");
            }
            return letters;
        }

        inline SqType guessSqType(const std::vector<std::string> &letters) {
            for (auto &type : {DNA, DNA_CLN, RNA, RNA_CLN, AMI, AMI_CLN}) {
                if (getStandardLetters(type) == letters) return type;
            }
            return UNT;
        }

        inline std::string getSqTypeAbbr(const SqType &type) {
            switch (type) {
                case AMI:       return "ami";
                case AMI_CLN:   return "ami_cln";
                case DNA:       return "dna";
                case DNA_CLN:   return "dna_cln";
                case RNA:       return "rna";
                case RNA_CLN:   return "rna_cln";
                case UNT:       return "unt";
                case ATP:       return "atp";
                case ENC:       return "enc";
                default:        throw std::exception();
            }
        }

        inline std::vector<std::string> getClassStringVector(const SqType &type) {
            switch (type) {
                case AMI:
                    return {"amisq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case AMI_CLN:
                    return {"amisq", "clnsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case DNA:
                    return {"dnasq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case DNA_CLN:
                    return {"dnasq", "clnsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case RNA:
                    return {"rnasq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case RNA_CLN:
                    return {"amisq", "clnsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case UNT:
                    return {"untsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case ATP:
                    return {"atpsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case ENC:
                    return {"encsq", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                default:
                    return {};
            }
        }

        inline SqType getSqType(const Rcpp::StringVector &classVector) {
            std::string type = getScalarStringValue(classVector);
            if (type == "amisq") {
                if (getScalarStringValue(classVector, 1) == "clnsq") return AMI_CLN;
                return AMI;
            } else if (type == "dnasq") {
                if (getScalarStringValue(classVector, 1) == "clnsq") return DNA_CLN;
                return DNA;
            } else if (type == "rnasq") {
                if (getScalarStringValue(classVector, 1) == "clnsq") return RNA_CLN;
                return RNA;
            } else if (type == "untsq") return UNT;
            else if (type == "atpsq") return ATP;
            else if (type == "encsq") return ENC;
            else throw std::invalid_argument("Object does not have a proper sq subtype!");
        }
    
        inline SqType get_sq_type_from_abbr(const Rcpp::StringVector &type_vector) {
            std::string type = getScalarStringValue(type_vector);
            if (type == "ami_bsc") return AMI_CLN;
            else if (type == "ami_ext") return AMI;
            else if (type == "dna_bsc") return DNA_CLN;
            else if (type == "dna_ext") return DNA;
            else if (type == "rna_bsc") return RNA_CLN;
            else if (type == "rna_ext") return RNA;
            else if (type == "unt") return UNT;
            else if (type == "atp") return ATP;
            else if (type == "enc") return ENC;
            else throw std::invalid_argument("404: type doesn't exist");
        }
    }
}

#endif //TIDYSQ_ALPHABET_H
