#ifndef TIDYSQ_ALPHABETSTD_H
#define TIDYSQ_ALPHABETSTD_H

#include <vector>
#include <string>
#include "general.h"

namespace tidysq {
    template<>
    class Alphabet<STD> : public std::vector<std::string> {
        using std::vector<std::string>::vector;
    };
}

#endif //TIDYSQ_ALPHABETSTD_H
