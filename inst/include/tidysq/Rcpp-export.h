#pragma once

#include <tidysq/internal/MotifFrame.h>
#include "tidysq/Sq.h"
#include "tidysq/ProtoSq.h"
#include "tidysq/find_motifs.h"
#include "tidysq/io/read_fasta.h"

namespace tidysq {
    inline Rcpp::StringVector export_to_R(const Alphabet &alphabet) {
        Rcpp::StringVector ret = util::convert_string_vector(util::convert_map_to_vector(alphabet.value_to_letter_));
        ret.attr("type") = util::sq_type_abbr_for_type(alphabet.type_);
        ret.attr("class") = Rcpp::StringVector{"sq_alphabet", "character", "vctrs_vctr"};
        return ret;
    }

    inline Rcpp::List export_to_R(const Sq<RCPP_IT> &sq) {
        Rcpp::List ret = sq.content_;
        ret.attr("alphabet") = export_to_R(sq.alphabet_);
        ret.attr("class") = util::sq_R_class_for_sq_type(sq.type());
        ret.attr("ptype") = Rcpp::RawVector{};
        return ret;
    }

    template<typename INTERNAL, typename PROTO>
    inline typename ProtoSq<INTERNAL, PROTO>::ContentStorageType export_to_R(const ProtoSq<INTERNAL, PROTO> &proto_sq) {
        return proto_sq.content_;
    }

    inline Rcpp::DataFrame export_to_R(const internal::NamedSqibble<RCPP_IT> &sqibble) {
        auto ret = Rcpp::DataFrame::create(Rcpp::Named("sq") = export_to_R(std::get<0>(sqibble)),
                                           Rcpp::Named("name") = util::convert_string_vector(std::get<1>(sqibble)));
        ret.attr("class") = Rcpp::StringVector{"tbl_df", "tbl", "data.frame"};
        return ret;
    }


    inline Rcpp::DataFrame export_to_R(const internal::MotifFrame<RCPP_IT> &found_motifs) {
        auto ret = Rcpp::DataFrame::create(
                Rcpp::Named("names", found_motifs.names()),
                Rcpp::Named("found", export_to_R(found_motifs.found())),
                Rcpp::Named("sought", found_motifs.sought()),
                Rcpp::Named("start", Rcpp::IntegerVector(Rcpp::wrap(found_motifs.start())) + 1),
                Rcpp::Named("end", Rcpp::IntegerVector(Rcpp::wrap(found_motifs.end())) + 1));
        ret.attr("class") = Rcpp::StringVector{"tbl_df", "tbl", "data.frame"};
        return ret;
    }
}