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

    Alphabet alphabet(std::vector<std::string>{"A", "C", "G", "T", "-"}, "!");

    SqProto<STD, RAWS> sq_proto(
        std::vector<SequenceProto<STD, RAWS>>{sequence_proto_1,
                                              sequence_proto_2,
                                              sequence_proto_3},
        alphabet);

    Sq<STD> sq = sq_proto.pack<STD>();

  }
  
  test_that("pack RCPP RAWS 1") {
    SequenceProto<RCPP, RAWS> sequence_proto_1{'A', 'C', 'T', 'G'};
    SequenceProto<RCPP, RAWS> sequence_proto_2{'A', 'A'};
    SequenceProto<RCPP, RAWS> sequence_proto_3{'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'};
    
    Rcpp::StringVector alph_content{"A", "C", "G", "T", "-"};
    alph_content.attr("na_character") = "!";
    
    Alphabet alphabet(alph_content);
    
    SqProto<RCPP, RAWS> sq_proto(
        Rcpp::List{sequence_proto_1,
                   sequence_proto_2,
                   sequence_proto_3},
        alphabet);

    Sq<RCPP> sq = sq_proto.pack<RCPP>();
  }
}