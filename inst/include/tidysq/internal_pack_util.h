#ifndef TIDYSQ_INTERNAL_PACK_UTIL_H
#define TIDYSQ_INTERNAL_PACK_UTIL_H

#include <cmath>
#include "types_general.h"

namespace tidysq {
    namespace internal {
        template<typename ALPHABET>
        constexpr inline sizealph getAlphabetSize(const ALPHABET& alphabet) {
            return ceil(log2((double) alphabet.size() + 2));
        }

        template<typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        constexpr inline lensq getPackedLength(const SEQUENCE_UNPACKED& unpacked, const ALPHABET& alphabet) {
            return (getAlphabetSize(alphabet) * unpacked.size() + 7) / 8;
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED reserveSpaceForPacked(const SEQUENCE_UNPACKED& unpacked, const ALPHABET& alphabet) {
            return SEQUENCE_PACKED(getPackedLength(unpacked, alphabet), unpacked.size());
        }
    }
}

#endif //TIDYSQ_INTERNAL_PACK_UTIL_H
