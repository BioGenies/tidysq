#pragma once

#include "tidysq/types/Sq.h"
#include "tidysq/types/ProtoSq.h"

namespace tidysq {
    inline Rcpp::StringVector export_to_R(const Alphabet &alphabet) {
        Rcpp::StringVector ret = util::convertStringVector(alphabet.letters_);
        ret.attr("type") = util::sq_type_abbr_for_type(alphabet.type_);
        ret.attr("class") = Rcpp::StringVector{"sq_alphabet", "character", "vctrs_vctr"};
        return ret;
    }

    inline Rcpp::List export_to_R(const Sq<RCPP> &sq) {
        Rcpp::List ret = sq.content_;
        ret.attr("alphabet") = export_to_R(sq.alphabet_);
        ret.attr("class") = util::sq_R_class_for_sq_type(sq.type());
        ret.attr("ptype") = Rcpp::RawVector{};
        return ret;
    }

    template<InternalType INTERNAL, ProtoType PROTO>
    inline typename ProtoSq<INTERNAL, PROTO>::ContentType export_to_R(const ProtoSq<INTERNAL, PROTO> &proto_sq) {
        return proto_sq.content_;
    }
}