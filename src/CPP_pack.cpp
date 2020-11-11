#include "tidysq/tidysq-includes.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(const Rcpp::List &proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<RAWS>(proto, alphabet, NA_letter)
                    .pack<RCPP>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(const Rcpp::List& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<INTS>(proto, alphabet, NA_letter)
                    .pack<RCPP>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(const Rcpp::List& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<STRINGS>(proto, alphabet, NA_letter)
                    .pack<RCPP>());
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(const Rcpp::StringVector& proto, const Rcpp::StringVector& alphabet, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_proto_from_R<STRING>(proto, alphabet, NA_letter)
                    .pack<RCPP>());
}