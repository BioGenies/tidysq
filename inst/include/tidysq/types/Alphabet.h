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
#include "tidysq/util/common.h"

namespace tidysq {
    namespace util {
        SqType guess_sq_type(const std::vector<Letter> &letters);

        Letter default_NA_letter();

        std::vector<Letter> standard_letters_for_type(const SqType &type);

        std::vector<std::string> sq_R_style_class_for_type(const SqType &type);

        SqType sq_type_from_R_style_class(const Rcpp::StringVector &classVector);

        std::string sq_type_abbr_for_type(const SqType &type);

        SqType sq_type_for_abbr(const Rcpp::StringVector &type_vector);
    }

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

        Rcpp::StringVector export_letters() {
            Rcpp::StringVector ret(letters_.size());
            auto iterator_in = letters_.begin();
            auto iterator_out = ret.begin();
            while (iterator_in != letters_.end()) {
                *iterator_out = *iterator_in;
                iterator_in++;
                iterator_out++;
            }
            return ret;
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
                      Alphabet(util::standard_letters_for_type(type),
                              type,
                              NA_letter) {};

        Alphabet(const Rcpp::StringVector &letters,
                 const SqType &type,
                 const Rcpp::StringVector &NA_letter = util::default_NA_letter()) :
                Alphabet(util::convertStringVector(letters),
                         type,
                         util::getScalarStringValue(NA_letter)) {};

        Alphabet(const Rcpp::List::const_AttributeProxy &letters,
                 const SqType &type,
                 const Rcpp::StringVector &NA_letter = util::default_NA_letter()) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(letters),
                         type,
                         util::getScalarStringValue(NA_letter)) {};

        explicit Alphabet(std::vector<Letter> letters,
                          const Letter &NA_letter = util::default_NA_letter()) :
                Alphabet(letters,
                         util::guess_sq_type(letters),
                         NA_letter) {};

        explicit Alphabet(const Rcpp::StringVector &letters,
                          const Rcpp::StringVector &NA_letter = util::default_NA_letter()) :
                Alphabet(letters,
                         util::sq_type_for_abbr(letters.attr("type")),
                         NA_letter) {};

        Alphabet(const Alphabet &other) = default;

        Alphabet(Alphabet &&other) noexcept = default;

        Rcpp::StringVector export_to_R() {
            Rcpp::StringVector ret(export_letters());
            ret.attr("type") = util::sq_type_abbr_for_type(type_);
            ret.attr("class") = Rcpp::StringVector{"sq_alphabet", "character", "vctrs_vctr"};
            return ret;
        }

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
    };

    namespace util {
        inline Letter default_NA_letter() {
            return "!";
        }

        inline std::vector<Letter> standard_letters_for_type(const SqType &type) {
            std::vector<Letter> letters;
            switch (type) {
                case AMI_EXT:
                    letters = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                               "S", "T", "U", "V", "W", "X", "Y", "Z", "-", "*"};
                    break;
                case AMI_BSC:
                    letters = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                               "W", "Y", "-", "*"};
                    break;
                case DNA_EXT:
                    letters = {"A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                    break;
                case DNA_BSC:
                    letters = {"A", "C", "G", "T", "-"};
                    break;
                case RNA_EXT:
                    letters = {"A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                    break;
                case RNA_BSC:
                    letters = {"A", "C", "G", "U", "-"};
                    break;
                default:
                    throw std::invalid_argument("Provided type does not have a predefined standard alphabet!");
            }
            return letters;
        }

        inline SqType guess_sq_type(const std::vector<Letter> &letters) {
            for (auto &type : {DNA_EXT, DNA_BSC, RNA_EXT, RNA_BSC, AMI_EXT, AMI_BSC}) {
                if (standard_letters_for_type(type) == letters) return type;
            }
            return UNT;
        }

        inline std::string sq_type_abbr_for_type(const SqType &type) {
            switch (type) {
                case AMI_EXT:   return "ami_ext";
                case AMI_BSC:   return "ami_bsc";
                case DNA_EXT:   return "dna_ext";
                case DNA_BSC:   return "dna_bsc";
                case RNA_EXT:   return "rna_ext";
                case RNA_BSC:   return "rna_bsc";
                case UNT:       return "unt";
                case ATP:       return "atp";
                case ENC:       return "enc";
                default:        throw std::exception();
            }
        }

        inline std::vector<std::string> sq_R_style_class_for_type(const SqType &type) {
            switch (type) {
                case AMI_EXT:
                    return {"sq_ami_ext", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case AMI_BSC:
                    return {"sq_ami_bsc", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case DNA_EXT:
                    return {"sq_dna_ext", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case DNA_BSC:
                    return {"sq_dna_bsc", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case RNA_EXT:
                    return {"sq_rna_ext", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case RNA_BSC:
                    return {"sq_ami_bsc", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case UNT:
                    return {"sq_unt", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case ATP:
                    return {"sq_atp", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                case ENC:
                    return {"sq_enc", "sq", "vctrs_list_of", "vctrs_vctr", "list"};
                default:
                    return {};
            }
        }

        inline SqType sq_type_from_R_style_class(const Rcpp::StringVector &classVector) {
            std::string type = getScalarStringValue(classVector);
            if (type == "sq_ami_bsc") return AMI_BSC;
            if (type == "sq_ami_ext") return AMI_EXT;
            if (type == "sq_dna_bsc") return DNA_BSC;
            if (type == "sq_dna_ext") return DNA_EXT;
            if (type == "sq_rna_bsc") return RNA_BSC;
            if (type == "sq_rna_ext") return RNA_EXT;
            if (type == "sq_unt") return UNT;
            if (type == "sq_atp") return ATP;
            if (type == "sq_enc") return ENC;
            throw std::invalid_argument("Object does not have a proper sq subtype!");
        }
    
        inline SqType sq_type_for_abbr(const Rcpp::StringVector &type_vector) {
            std::string type = getScalarStringValue(type_vector);
            if (type == "ami_bsc") return AMI_BSC;
            if (type == "ami_ext") return AMI_EXT;
            if (type == "dna_bsc") return DNA_BSC;
            if (type == "dna_ext") return DNA_EXT;
            if (type == "rna_bsc") return RNA_BSC;
            if (type == "rna_ext") return RNA_EXT;
            if (type == "unt") return UNT;
            if (type == "atp") return ATP;
            if (type == "enc") return ENC;
            else throw std::invalid_argument("404: type doesn't exist");
        }
    }
}

#endif //TIDYSQ_ALPHABET_H