#ifndef TIDYSQ_GENERAL_H
#define TIDYSQ_GENERAL_H

#include <string>

namespace tidysq {
    typedef unsigned long long int  LenSq;
    typedef unsigned char           ElemPacked;
    typedef unsigned char           ElemRaws;
    typedef unsigned short int      ElemInts;
    typedef std::string             ElemStrings;
    typedef char                    ElemString;
    typedef unsigned short int      AlphSize;
    typedef unsigned short int      LetValue;

    enum InternalType {
        STD,
        RCPP
    };

    enum ProtoType {
        RAWS,
        INTS,
        STRINGS,
        STRING
    };

    enum SqType {
        AMI,
        AMI_CLN,
        DNA,
        DNA_CLN,
        RNA,
        RNA_CLN,
        UNT,
        ATP,
        ENC
    };
}

#endif //TIDYSQ_GENERAL_H
