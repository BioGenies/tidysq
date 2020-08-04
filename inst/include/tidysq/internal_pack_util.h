#ifndef TIDYSQ_INTERNAL_PACK_UTIL_H
#define TIDYSQ_INTERNAL_PACK_UTIL_H

#include <cmath>
#include "types_Sequence.h"
#include "types_RcppSequence.h"


namespace tidysq {
    namespace internal {
        template<typename ALPHABET>
        inline sizealph getAlphabetSize(const ALPHABET& alphabet) {
            return ceil(log2((double) alphabet.size() + 2));
        }

        template<typename SEQUENCE_UNPACKED,
                typename ALPHABET>
        inline lensq getPackedLength(const SEQUENCE_UNPACKED& unpacked, const ALPHABET& alphabet) {
            return (getAlphabetSize(alphabet) * unpacked.size() + 7) / 8;
        }

        template<typename SEQUENCE_PACKED>
        inline lensq getOriginalLength(const SEQUENCE_PACKED& packed);

        template<>
        inline lensq getOriginalLength<Sequence>(const Sequence& packed) {
            return packed.originalLength();
        }

        template<>
        inline lensq getOriginalLength<RcppSequence>(const RcppSequence& packed) {
            return packed.attr("original_length");
        }

        template<typename SEQUENCE_UNPACKED,
                typename SEQUENCE_PACKED,
                typename ALPHABET>
        inline SEQUENCE_PACKED reserveSpaceForPacked(const SEQUENCE_UNPACKED& unpacked, const ALPHABET& alphabet) {
            return SEQUENCE_PACKED(getPackedLength(unpacked, alphabet), unpacked.size());
        }

        template<typename SEQUENCE_PACKED,
                typename SEQUENCE_UNPACKED>
        inline SEQUENCE_UNPACKED reserveSpaceForUnpacked(const SEQUENCE_PACKED& packed) {
            return SEQUENCE_UNPACKED(getOriginalLength(packed));
        }
    }
}

#endif //TIDYSQ_INTERNAL_PACK_UTIL_H
