#ifndef TIDYSQ_SEQUENCESTD_H
#define TIDYSQ_SEQUENCESTD_H

#include <vector>
#include <string>

#include "general.h"

namespace tidysq {
    template<>
    class Sequence<STD> : public std::pair<std::vector<unsigned char>, lensq> {
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
}

#endif //TIDYSQ_SEQUENCESTD_H
