#pragma once

#include <map>
#include <algorithm>

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    const internal::AmbiguousDict AMBIGUOUS_AMINO_MAP = {
            {"B", {"B", "D", "N"}},
            {"J", {"J", "I", "L"}},
            {"Z", {"Z", "E", "Q"}},
            {"X", {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                          "S", "T", "U", "V", "W", "X", "Y", "Z"}}
    };

    const internal::AmbiguousDict AMBIGUOUS_DNA_MAP = {
            {"W", {"W", "A", "T"}},
            {"S", {"S", "C", "G"}},
            {"M", {"M", "A", "C"}},
            {"K", {"K", "G", "T"}},
            {"R", {"R", "A", "G"}},
            {"Y", {"Y", "C", "T"}},
            {"B", {"B", "S", "K", "Y", "C", "G", "T"}},
            {"D", {"D", "W", "K", "R", "A", "G", "T"}},
            {"H", {"H", "W", "M", "Y", "A", "C", "T"}},
            {"V", {"V", "S", "M", "R", "A", "C", "G"}},
            {"N", {"A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"}}
    };

    const internal::AmbiguousDict AMBIGUOUS_RNA_MAP = {
            {"W", {"W", "A", "U"}},
            {"S", {"S", "C", "G"}},
            {"M", {"M", "A", "C"}},
            {"K", {"K", "G", "U"}},
            {"R", {"R", "A", "G"}},
            {"Y", {"Y", "C", "U"}},
            {"B", {"B", "S", "K", "Y", "C", "G", "U"}},
            {"D", {"D", "W", "K", "R", "A", "G", "U"}},
            {"H", {"H", "W", "M", "Y", "A", "C", "U"}},
            {"V", {"V", "S", "M", "R", "A", "C", "G"}},
            {"N", {"A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"}}
    };
}
