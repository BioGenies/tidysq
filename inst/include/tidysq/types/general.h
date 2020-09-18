#ifndef TIDYSQ_GENERAL_H
#define TIDYSQ_GENERAL_H

namespace tidysq {
    typedef unsigned long long int lensq;
    typedef unsigned char ElemRaws;
    typedef unsigned short int ElemInts;
    typedef std::string ElemStrings;
    typedef char ElemString;
    typedef unsigned short int sizealph;
    typedef unsigned short int letvalue;

    enum InternalType {
        ANY_INTERNAL,
        STD,
        RCPP
    };

    enum ProtoType {
        ANY_PROTO,
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

    template<InternalType INTERNAL = STD>
    class Sequence;

    template<InternalType INTERNAL = STD>
    class Sq;

    template<InternalType INTERNAL = STD,
            ProtoType PROTO = RAWS>
    class SequenceProto;

    template<InternalType INTERNAL = STD,
            ProtoType PROTO = RAWS>
    class SqProto;

    class Alphabet;
}

#endif //TIDYSQ_GENERAL_H
