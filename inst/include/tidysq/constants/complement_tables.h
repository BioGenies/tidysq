#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    const ComplementTable BSC_COMPLEMENT_TABLE = {
            {0u, 3u}, {1u, 2u}, {2u, 1u}, {3u, 0u}
    };

    const ComplementTable EXT_COMPLEMENT_TABLE = {
            {0u, 3u}, {1u, 2u}, {2u, 1u}, {3u, 0u},
            {4u, 4u}, {5u, 5u}, {6u, 7u}, {7u, 6u}, {8u, 9u}, {9u, 8u},
            {10u, 13u}, {11u, 12u}, {12u, 11u}, {13u, 10u},
            {14u, 14u}
    };
}
