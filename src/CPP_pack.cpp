#include "tidysq/tidysq-includes.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(const Rcpp::List &proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<RAWS_PT>(proto, alphabet, NA_letter)
                    .pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(const Rcpp::List& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<INTS_PT>(proto, alphabet, NA_letter)
                    .pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(const Rcpp::List& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<STRINGS_PT>(proto, alphabet, NA_letter)
                    .pack<RCPP_IT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(const Rcpp::StringVector& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<STRING_PT>(proto, alphabet, NA_letter)
                    .pack<RCPP_IT>());
}