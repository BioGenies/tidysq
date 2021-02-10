#pragma once

#include <unordered_map>

#include "tidysq/sq-types.h"

namespace tidysq::constants {
    const Letter DEFAULT_NA_LETTER{"!"};
    const bool DEFAULT_IGNORE_CASE = false;

    const std::unordered_map<SqType, const std::vector<Letter>>  STANDARD_LETTERS{
            {DNA_BSC, {"A", "C", "G", "T", "-"}},
            {DNA_EXT, {"A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"}},
            {RNA_BSC, {"A", "C", "G", "U", "-"}},
            {RNA_EXT, {"A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"}},
            {AMI_BSC, {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                              "W", "Y", "-", "*"}},
            {AMI_EXT, {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                              "S", "T", "U", "V", "W", "X", "Y", "Z", "-", "*"}}
    };
}
