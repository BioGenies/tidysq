#include "tidysq/internal/pack.h"

#include "tidysq/internal/bit_macros_defs.h"

namespace tidysq::v2::internal {
    template<BitIndex LEN, typename ITER_CONST_IN, typename ITER_OUT>
    void pack_3(ITER_CONST_IN &it_in, ITER_OUT &it_out) {
        LetterVal v1, v2, v3, v4;
        FETCH(1)
        FETCH_OR_ZERO(2, 1)
        FETCH_OR_ZERO(3, 2)
        ALIGN_3(1, 0, 2, 3, 3, 6)
        FETCH_OR_RETURN(1, 3)
        FETCH_OR_ZERO(2, 4)
        FETCH_OR_ZERO(4, 5)
        ALIGN_4(3, 2, 1, 1, 2, 4, 4, 7)
        FETCH_OR_RETURN(1, 6)
        FETCH_OR_ZERO(2, 7)
        ALIGN_3(4, 1, 1, 2, 2, 5)
    }

    template<BitIndex LEN, typename ITER_CONST_IN, typename ITER_OUT>
    void pack_5(ITER_CONST_IN &it_in, ITER_OUT &it_out) {
        LetterVal v1, v2, v3;
        FETCH(1)
        FETCH_OR_ZERO(2, 1)
        ALIGN_2(1, 0, 2, 5)
        FETCH_OR_RETURN(1, 2)
        FETCH_OR_ZERO(3, 3)
        ALIGN_3(2, 3, 1, 2, 3, 7)
        FETCH_OR_RETURN(1, 4)
        ALIGN_2(3, 1, 1, 4)
        FETCH_OR_RETURN(2, 5)
        FETCH_OR_ZERO(3, 6)
        ALIGN_3(1, 4, 2, 1, 3, 6)
        FETCH_OR_RETURN(1, 7)
        ALIGN_2(3, 2, 1, 3)
    }
}

#include "tidysq/internal/bit_macros_undefs.h"

template void tidysq::v2::internal::pack_3<8, std::vector<unsigned char>::const_iterator , std::vector<unsigned char>::iterator>(
        std::vector<unsigned char>::const_iterator &it_in, std::vector<unsigned char>::iterator &it_out);

template void tidysq::v2::internal::pack_5<8, std::vector<unsigned char>::const_iterator , std::vector<unsigned char>::iterator>(
        std::vector<unsigned char>::const_iterator &it_in, std::vector<unsigned char>::iterator &it_out);
