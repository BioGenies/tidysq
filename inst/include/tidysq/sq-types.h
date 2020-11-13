#pragma once

#include <vector>
#include <stdexcept>

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/constants/standard_letters.h"

namespace tidysq {
    typedef std::string SqTypeAbbr;
    typedef std::string RClass;

    namespace util {
        inline bool has_standard_alphabet(const SqType &type) {
            return std::set{AMI_EXT, AMI_BSC, DNA_EXT, DNA_BSC, RNA_EXT, RNA_BSC}.count(type);
        }

        inline std::vector<Letter> standard_letters_for_sq_type(const SqType &type) {
            try {
                return constants::STANDARD_LETTERS.at(type);
            } catch (const std::out_of_range &e) {
                throw std::invalid_argument("Provided R_class does not have a predefined standard alphabet!");
            }
        }

        inline SqType guess_sq_type_from_letters(const std::vector<Letter> &letters) {
            for (auto &type : {DNA_BSC, RNA_BSC, AMI_BSC, DNA_EXT, RNA_EXT, AMI_EXT}) {
                auto standard_alph = standard_letters_for_sq_type(type);
                if (std::all_of(letters.begin(), letters.end(), [=](const Letter &letter) {
                    return std::any_of(standard_alph.begin(), standard_alph.end(), [=](const Letter &other) {
                        return other == letter;
                    });
                })) {
                    return type;
                }
            }
            return UNT;
        }

        inline SqTypeAbbr sq_type_abbr_for_type(const SqType &type) {
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
                default:        throw std::invalid_argument("Provided R_class does not exist!");
            }
        }

        inline std::vector<RClass> sq_R_class_for_sq_type(const SqType &type) {
            return {"sq_" + sq_type_abbr_for_type(type),  "sq", "vctrs_list_of", "vctrs_vctr", "list"};
        }

        inline SqType sq_type_for_R_class(const RClass &R_class) {
            if (R_class == "sq_ami_bsc") return AMI_BSC;
            if (R_class == "sq_ami_ext") return AMI_EXT;
            if (R_class == "sq_dna_bsc") return DNA_BSC;
            if (R_class == "sq_dna_ext") return DNA_EXT;
            if (R_class == "sq_rna_bsc") return RNA_BSC;
            if (R_class == "sq_rna_ext") return RNA_EXT;
            if (R_class == "sq_unt") return UNT;
            if (R_class == "sq_atp") return ATP;
            if (R_class == "sq_enc") return ENC;
            throw std::invalid_argument("Object does not have a proper sq subtype!");
        }

        inline SqType sq_type_for_sq_type_abbr(const SqTypeAbbr &type_abbr) {
            if (type_abbr == "ami_bsc") return AMI_BSC;
            if (type_abbr == "ami_ext") return AMI_EXT;
            if (type_abbr == "dna_bsc") return DNA_BSC;
            if (type_abbr == "dna_ext") return DNA_EXT;
            if (type_abbr == "rna_bsc") return RNA_BSC;
            if (type_abbr == "rna_ext") return RNA_EXT;
            if (type_abbr == "unt") return UNT;
            if (type_abbr == "atp") return ATP;
            if (type_abbr == "enc") return ENC;
            else throw std::invalid_argument("404: R_class doesn't exist");
        }
    }
}