#ifndef TIDYSQ_COMMON_H
#define TIDYSQ_COMMON_H

#include <vector>
#include <string>
#include <Rcpp.h>


namespace tidysq::util {
    inline std::vector<std::string> convertStringVector(const Rcpp::StringVector &vector) {
        std::vector<std::string> ret(vector.size());
        auto iterIn = vector.begin();
        auto iterOut = ret.begin();
        while (iterIn != vector.end()) {
            *iterOut = *iterIn;
            iterIn++;
            iterOut++;
        }
        return ret;
    }

    inline Rcpp::StringVector convertStringVector(const std::vector<std::string> &vector) {
        Rcpp::StringVector ret(vector.size());
        auto iterIn = vector.begin();
        auto iterOut = ret.begin();
        while (iterIn != vector.end()) {
            *iterOut = *iterIn;
            iterIn++;
            iterOut++;
        }
        return ret;
    }

    inline std::string getScalarStringValue(const Rcpp::StringVector &vector, const unsigned int index = 0) {
        return Rcpp::as<std::string>(vector[index]);
    }
}

#endif //TIDYSQ_COMMON_H
