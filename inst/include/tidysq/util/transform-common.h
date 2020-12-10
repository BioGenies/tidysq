#pragma once

#include <vector>
#include <string>
#include <Rcpp.h>

namespace tidysq::util {
    inline std::vector<std::string> convert_string_vector(const Rcpp::StringVector &vector) {
        std::vector<std::string> ret(vector.size());
        std::copy(vector.begin(), vector.end(), ret.begin());
        return ret;
    }

    inline Rcpp::StringVector convert_string_vector(const std::vector<std::string> &vector) {
        Rcpp::StringVector ret(vector.size());
        std::copy(vector.begin(), vector.end(), ret.begin());
        return ret;
    }

    inline std::string convert_to_scalar(const Rcpp::StringVector &vector, const unsigned int index = 0) {
        return Rcpp::as<std::string>(vector[index]);
    }

    inline int convert_to_scalar(const Rcpp::IntegerVector &vector, const unsigned int index = 0) {
        return vector[index];
    }

    inline bool convert_to_scalar(const Rcpp::LogicalVector &vector, const unsigned int index = 0) {
        return vector[index];
    }

    template<typename T>
    inline std::vector<T> convert_set_to_vector(const std::set<T> &set) {
        std::vector<T> ret(set.size());
        std::copy(set.begin(), set.end(), ret.begin());
        return ret;
    }

    inline std::vector<Letter> convert_map_to_vector(const std::unordered_map<LetterValue, const Letter> &value_to_letter) {
        std::vector<Letter> ret(value_to_letter.size());
        for (unsigned short i = 0; i < value_to_letter.size(); i++) {
            ret[i] = value_to_letter.at(i);
        }
        return ret;
    }

    inline LenSq convert_sample_size(const Rcpp::NumericVector &sample_size) {
        return Rcpp::traits::is_infinite<REALSXP>(sample_size[0]) ? R_XLEN_T_MAX : sample_size[0];
    }
}
