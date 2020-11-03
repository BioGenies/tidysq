#ifndef TIDYSQ_GENERAL_H
#define TIDYSQ_GENERAL_H

#include <string>

namespace tidysq {
    typedef unsigned long long int  LenSq;
    typedef unsigned char           ElementPacked;
    typedef unsigned char           ElementRaws;
    typedef unsigned short int      ElementInts;
    typedef std::string             ElementStrings;
    typedef char                    ElementStringSimple;
    typedef std::string             ElementStringMultichar;
    typedef unsigned short int      AlphSize;
    typedef unsigned short int      LetterValue;
    typedef std::string             Letter;
    typedef char                    SimpleLetter;

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
}

#endif //TIDYSQ_GENERAL_H
