#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    inline Alphabet import_alphabet_from_R(const Rcpp::StringVector &letters, const Rcpp::StringVector &NA_letter) {
        return Alphabet(util::convert_string_vector(letters),
                        util::sq_type_for_sq_type_abbr(letters.attr("type")),
                        util::get_scalar_string_value(NA_letter));
    }

    inline Sq<RCPP_IT> import_from_R(const Rcpp::List &sq, const Rcpp::StringVector &NA_letter) {
        if (!sq.hasAttribute("alphabet")) throw std::invalid_argument("Sq object should have 'alphabet' attribute.");
        Rcpp::StringVector alphabet = sq.attr("alphabet");
        return Sq<RCPP_IT>(sq, import_alphabet_from_R(alphabet, NA_letter));
    }

    template<typename PROTO>
    ProtoSq<RCPP_IT, PROTO> import_proto_from_R(const typename ProtoSq<RCPP_IT, PROTO>::ContentType &proto,
                                             const Rcpp::StringVector &alphabet,
                                             const Rcpp::StringVector &NA_letter) {
        return ProtoSq<RCPP_IT, PROTO>(proto, import_alphabet_from_R(alphabet, NA_letter));
    }
}