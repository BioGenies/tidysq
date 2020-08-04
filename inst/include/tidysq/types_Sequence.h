#ifndef TIDYSQ_TYPES_SEQUENCE_H
#define TIDYSQ_TYPES_SEQUENCE_H

#include <vector>
#include <string>

namespace tidysq {
    typedef std::vector<unsigned char> Sequence;

    typedef std::vector<unsigned char> SequenceProtoRaw;
    typedef std::vector<int> SequenceProtoInteger;
    typedef std::string SequenceProtoString;
    typedef std::vector<std::string> SequenceProtoStrings;

    typedef std::vector<std::string> Alphabet;
}

#endif //TIDYSQ_TYPES_SEQUENCE_H
