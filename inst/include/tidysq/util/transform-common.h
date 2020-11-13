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

    inline std::string get_scalar_string_value(const Rcpp::StringVector &vector, const unsigned int index = 0) {
        return Rcpp::as<std::string>(vector[index]);
    }

    template<typename T>
    inline std::vector<T> convert_set_to_vector(const std::set<T> &set) {
        std::vector<T> ret(set.size());
        std::copy(set.begin(), set.end(), ret.begin());
        return ret;
    }
}
