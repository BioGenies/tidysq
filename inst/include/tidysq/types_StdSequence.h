#ifndef TIDYSQ_TYPES_STDSEQUENCE_H
#define TIDYSQ_TYPES_STDSEQUENCE_H

#include <vector>
#include <string>

#include "types_general.h"

namespace tidysq {
    class StdSequence : public std::pair<std::vector<unsigned char>, lensq> {
        typedef std::pair<std::vector<unsigned char>, lensq> BaseType;
    public:
        typedef unsigned char ElementType;
        StdSequence() :
                StdSequence(0, 0) {};

        StdSequence(lensq packed_length, lensq original_length) :
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

    typedef std::vector<unsigned char> StdSequenceProtoRaw;
    typedef std::vector<int> StdSequenceProtoInteger;
    typedef std::string StdSequenceProtoString;
    typedef std::vector<std::string> StdSequenceProtoStrings;

    typedef std::vector<std::string> StdAlphabet;
}

#endif //TIDYSQ_TYPES_STDSEQUENCE_H
