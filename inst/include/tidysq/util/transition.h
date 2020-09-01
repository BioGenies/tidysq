#ifndef TIDYSQ_TRANSITION_H
#define TIDYSQ_TRANSITION_H

#include <string>
#include <Rcpp.h>

#include "tidysq/types/general.h"

namespace tidysq::util {
    template<int DUMMY>
    std::string getNACharacterAsString(const Rcpp::StringVector &alphabet) {
        return Rcpp::as<std::string>(Rcpp::as<Rcpp::StringVector>(alphabet.attr("na_letter"))[0]);
    }

    template<int DUMMY>
    std::vector<std::string> getClassStringVector(const SqType &type) {
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

    template<int DUMMY>
    SqType getSqType(const Rcpp::StringVector &classVector) {
        std::string type = getScalarStringValue<0>(classVector);
        if (type == "amisq") {
            if (getScalarStringValue<0>(classVector, 1) == "clnsq") return AMI_CLN;
            return AMI;
        } else if (type == "dnasq") {
            if (getScalarStringValue<0>(classVector, 1) == "clnsq") return DNA_CLN;
            return DNA;
        } else if (type == "rnasq") {
            if (getScalarStringValue<0>(classVector, 1) == "clnsq") return RNA_CLN;
            return RNA;
        } else if (type == "untsq") {
            return UNT;
        } else if (type == "atpsq") {
            return ATP;
        }
        return ENC;
    }
}

#endif //TIDYSQ_TRANSITION_H
