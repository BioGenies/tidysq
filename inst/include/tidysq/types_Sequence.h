#ifndef TIDYSQ_TYPES_SEQUENCE_H
#define TIDYSQ_TYPES_SEQUENCE_H

#include <vector>
#include <string>

#include "types_general.h"

namespace tidysq {
    class Sequence : public std::pair<std::vector<unsigned char>, lensq> {
        typedef std::pair<std::vector<unsigned char>, lensq> BaseType;
    public:
        typedef unsigned char ElementType;
        Sequence() :
                Sequence(0, 0) {};

        Sequence(lensq packed_length, lensq original_length) :
                BaseType(std::vector<ElementType>(packed_length), original_length) {};

        const ElementType &operator[](lensq index) const {
            return first[index];
        }

        ElementType &operator[](lensq index) {
            return first[index];
        }

        [[nodiscard]] const lensq &originalLength() const {
            return second;
        }

        lensq &originalLength() {
            return second;
        }

        [[nodiscard]] lensq size() const {
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
