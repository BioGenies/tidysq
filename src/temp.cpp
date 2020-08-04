#include <Rcpp.h>
#include "tidysq/types_RcppSq.h"
#include "tidysq/types_RcppSqProto.h"

using namespace Rcpp;
using namespace tidysq;

//' @export
//[[Rcpp::export]]
List tmpPack(List raws, StringVector alphabet) {
  RcppSqProto<RcppSequenceProtoRaw> sqProto(raws, alphabet);
  return sqProto.pack<RcppSq>().exportToR();
}
