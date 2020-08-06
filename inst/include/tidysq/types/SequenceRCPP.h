#ifndef TIDYSQ_SEQUENCERCPP_H
#define TIDYSQ_SEQUENCERCPP_H

#include <Rcpp.h>

#include "general.h"

namespace tidysq {
    template<>
    class Sequence<RCPP> : public Rcpp::RawVector {
        typedef Rcpp::RawVector BaseType;
    public:
        Sequence(lensq packed_length, lensq original_length) :
                BaseType(packed_length) {
            BaseType::attr("original_length") = original_length;
        }

        explicit Sequence(const Rcpp::RawVector& content) :
                BaseType(content) {
            if (!content.hasAttribute("original_length"))
                throw std::invalid_argument(R"("content" argument in Sequence should have "original_length" attribute!)");
        };

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
}

#endif //TIDYSQ_SEQUENCERCPP_H
