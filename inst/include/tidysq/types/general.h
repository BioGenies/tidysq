#ifndef TIDYSQ_GENERAL_H
#define TIDYSQ_GENERAL_H

namespace tidysq {
    typedef unsigned long long int lensq;
    typedef short sizealph;
    typedef short letvalue;

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

    template<InternalType INTERNAL = STD>
    class Sequence;

    template<InternalType INTERNAL = STD>
    class Sq;

    template<InternalType INTERNAL = STD>
    class Alphabet;

    template<InternalType INTERNAL = STD,
            ProtoType PROTO = RAWS>
    class SequenceProto;

    template<InternalType INTERNAL = STD,
            ProtoType PROTO = RAWS>
    class SqProto;
}

#endif //TIDYSQ_GENERAL_H
