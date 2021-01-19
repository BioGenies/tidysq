#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(const Rcpp::List &proto,
                         const Rcpp::StringVector& alphabet,
                         const tidysq::Letter& NA_letter,
                         const bool& ignore_case) {
    return export_to_R(import_proto_sq_from_R<RAWS_PT>(proto, alphabet, NA_letter, ignore_case).pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(const Rcpp::List& proto,
                         const Rcpp::StringVector& alphabet,
                         const tidysq::Letter& NA_letter,
                         const bool& ignore_case) {
    return export_to_R(import_proto_sq_from_R<INTS_PT>(proto, alphabet, NA_letter, ignore_case).pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(const Rcpp::List& proto,
                            const Rcpp::StringVector& alphabet,
                            const tidysq::Letter& NA_letter,
                            const bool& ignore_case) {
    return export_to_R(import_proto_sq_from_R<STRINGS_PT>(proto, alphabet, NA_letter, ignore_case).pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(const Rcpp::StringVector& proto,
                           const Rcpp::StringVector& alphabet,
                           const tidysq::Letter& NA_letter,
                           const bool& ignore_case) {
    return export_to_R(import_proto_sq_from_R<STRING_PT>(proto, alphabet, NA_letter, ignore_case).pack<RCPP_IT>());
}