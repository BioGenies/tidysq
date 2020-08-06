#ifndef TIDYSQ_SEQUENCEPROTOSTD_H
#define TIDYSQ_SEQUENCEPROTOSTD_H

#include <vector>
#include <string>
#include "general.h"

namespace tidysq {
    template<>
    class SequenceProto<STD, RAWS> : public std::vector<unsigned char> {
        using std::vector<unsigned char>::vector;
    };

    template<>
    class SequenceProto<STD, INTS> : public std::vector<short int> {
        using std::vector<short int>::vector;
    };

    template<>
    class SequenceProto<STD, STRINGS> : public std::vector<std::string> {
        using std::vector<std::string>::vector;
    };

    template<>
    class SequenceProto<STD, STRING> : public std::string {
        using std::string::string;
    };

}

#endif //TIDYSQ_SEQUENCEPROTOSTD_H
