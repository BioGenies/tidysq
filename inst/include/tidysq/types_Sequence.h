#ifndef TIDYSQ_TYPES_SEQUENCE_H
#define TIDYSQ_TYPES_SEQUENCE_H

#include <vector>
#include <string>

#include "types_general.h"

namespace tidysq {
    class Sequence : public std::pair<std::vector<unsigned char>, lensq> {
        typedef std::pair<std::vector<unsigned char>, lensq> BaseType;
    public:
        Sequence(lensq packed_length, lensq original_length) :
                BaseType(std::vector<unsigned char>(packed_length), original_length) {};

        const unsigned char &operator[](lensq index) const {
            return first[index];
        }

        unsigned char &operator[](lensq index) {
            return first[index];
        }

        const lensq &originalLength() const {
            return second;
        }

        lensq &originalLength() {
            return second;
        }

        lensq size() const {
            return first.size();
        }
    };

    typedef std::vector<unsigned char> SequenceProtoRaw;
    typedef std::vector<int> SequenceProtoInteger;
    typedef std::string SequenceProtoString;
    typedef std::vector<std::string> SequenceProtoStrings;

    typedef std::vector<std::string> Alphabet;
}

#endif //TIDYSQ_TYPES_SEQUENCE_H
