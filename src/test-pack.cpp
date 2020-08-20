#include <Rcpp.h>
#include <testthat.h>

#include "tidysq/types/general.h"
#include "tidysq/types/SqSTD.h"
#include "tidysq/types/SqProtoSTD.h"
#include "tidysq/types/SequenceSTD.h"

using namespace tidysq;

context("test_packing") {
  test_that("pack STD RAWS 1") {
    SequenceProto<STD, RAWS> sequence_proto_1{'A', 'C', 'T', 'G'};
    SequenceProto<STD, RAWS> sequence_proto_2{'A', 'A'};
    SequenceProto<STD, RAWS> sequence_proto_3{'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'};
    
    Alphabet<STD> alphabet(std::vector<std::string>{"A", "C", "G", "T", "-"}, "!");
    
    SqProto<STD, RAWS> sq_proto(
        std::vector<SequenceProto<STD, RAWS>>{sequence_proto_1,
                                              sequence_proto_2,
                                              sequence_proto_3},
        alphabet);
    
    Sq<STD> sq = sq_proto.pack<STD>();

  }
}