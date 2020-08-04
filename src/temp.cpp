#include <Rcpp.h>
#include "tidysq/types_RcppSq.h"
#include "tidysq/types_RcppSqProto.h"

using namespace Rcpp;
using namespace tidysq;

//' @export
//[[Rcpp::export]]
List tmpPack(List raws, StringVector alphabet) {
  return RcppSqProto<RcppSequenceProtoRaw>(raws, alphabet)
    .pack<RcppSq>()
    .exportToR();
}

//' @export
//[[Rcpp::export]]
List tmpUnpack(List raws) {
  return RcppSq(raws)
    .unpack<RcppSqProto<RcppSequenceProtoRaw>>()
    .exportToR();
}