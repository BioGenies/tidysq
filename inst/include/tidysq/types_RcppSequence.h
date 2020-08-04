#ifndef TIDYSQ_TYPES_RCPPSEQUENCE_H
#define TIDYSQ_TYPES_RCPPSEQUENCE_H

#include <Rcpp.h>

#include "types_general.h"

namespace tidysq {
    class RcppSequence : public Rcpp::RawVector {
        typedef Rcpp::RawVector BaseType;
    public:
        RcppSequence(lensq packed_length, lensq original_length) :
                BaseType(packed_length) {
            BaseType::attr("original_length") = original_length;
        }

        BaseType::const_Proxy operator[] (lensq index) const {
            return BaseType::operator[](index);
        }

        BaseType::Proxy operator[] (lensq index) {
            return BaseType::operator[](index);
        }

        BaseType::const_AttributeProxy attr(const std::string& name) const {
            return BaseType::attr(name);
        }

        lensq size() const {
            return BaseType::size();
        }
    };

    typedef Rcpp::RawVector RcppSequenceProtoRaw;
    typedef Rcpp::IntegerVector RcppSequenceProtoInteger;
    typedef Rcpp::StringVector RcppSequenceProtoString;
    typedef Rcpp::StringVector RcppSequenceProtoStrings;

    typedef Rcpp::StringVector RcppAlphabet;
}

#endif //TIDYSQ_TYPES_RCPPSEQUENCE_H
