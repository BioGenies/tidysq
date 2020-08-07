#ifndef TIDYSQ_ALPHABETRCPP_H
#define TIDYSQ_ALPHABETRCPP_H

#include <Rcpp.h>

#include "general.h"

namespace tidysq {
    template<>
    class Alphabet<RCPP> : public Rcpp::StringVector {
    private:
        typedef Rcpp::StringVector BaseType;
        letvalue NAValue_;
        sizealph alphabetSize_;
        Rcpp::StringVector NALetter_;
    public:
        explicit Alphabet(const Rcpp::StringVector& letters) :
                BaseType(letters) {
            if (!letters.hasAttribute("na_character"))
                throw std::invalid_argument(R"("letters" argument in Alphabet should have "na_character" attribute!)");
            NALetter_ = letters.attr("na_character");
            alphabetSize_ = ceil(log2((double) letters.size() + 2));
            NAValue_ = pow(2, alphabetSize_) - 1;
        }

        explicit Alphabet(const Rcpp::List::const_AttributeProxy &letters) :
                Alphabet(Rcpp::as<Rcpp::StringVector>(letters)) {};

        [[nodiscard]] const letvalue &NAValue() const {
            return NAValue_;
        };

        letvalue &NAValue() {
            return NAValue_;
        }

        [[nodiscard]] Rcpp::String NALetter() const {
            return NALetter_[0];
        }

        [[nodiscard]] const sizealph &alphabetSize() const {
            return alphabetSize_;
        }
    };
}

#endif //TIDYSQ_ALPHABETRCPP_H
