#include "tidysq/internal/unpack_common.h"
#include "tidysq/internal/bit_macros_defs.h"

namespace tidysq::v2::internal {
    template<BitIndex LEN, typename ITER_CONST_IN, typename ITER_OUT>
    void unpack_3(ITER_CONST_IN &it_in, ITER_OUT &it_out) {
        LetterVal v1, v2;
        BitShift full_byte_trim = 7;
        FETCH(1)
        ALIGN_AND_TRIM_1(1, 0)
        ALIGN_AND_TRIM_1_OR_RETURN(1, 3, 1)
        FETCH_OR_RETURN(2, 2)
        ALIGN_AND_TRIM_2(1, 6, 3, 2, 2)
        ALIGN_AND_TRIM_1_OR_RETURN(2, 1, 3)
        ALIGN_AND_TRIM_1_OR_RETURN(2, 4, 4)
        FETCH_OR_RETURN(1, 5)
        ALIGN_AND_TRIM_2(2, 7, 1, 1, 1)
        ALIGN_AND_TRIM_1_OR_RETURN(1, 2, 6)
        ALIGN_AND_TRIM_1_OR_RETURN(1, 5, 7)
    }

    template<BitIndex LEN, typename ITER_CONST_IN, typename ITER_OUT>
    void unpack_5(ITER_CONST_IN &it_in, ITER_OUT &it_out) {
        LetterVal v1, v2;
        BitShift full_byte_trim = 31;
        FETCH(1)
        ALIGN_AND_TRIM_1(1, 0)
        FETCH_OR_RETURN(2, 1)
        ALIGN_AND_TRIM_2(1, 5, 7, 2, 3)
        ALIGN_AND_TRIM_1_OR_RETURN(2, 2, 2)
        FETCH_OR_RETURN(1, 3)
        ALIGN_AND_TRIM_2(2, 7, 1, 1, 1)
        FETCH_OR_RETURN(2, 4)
        ALIGN_AND_TRIM_2(1, 4, 15, 2, 4)
        ALIGN_AND_TRIM_1_OR_RETURN(2, 1, 5)
        FETCH_OR_RETURN(1, 6)
        ALIGN_AND_TRIM_2(2, 6, 3, 1, 2)
        ALIGN_AND_TRIM_1_OR_RETURN(1, 3, 7)
    }
}

#include "tidysq/internal/bit_macros_undefs.h"

template
void tidysq::v2::internal::unpack_3<8, std::vector<unsigned char>::const_iterator , std::vector<unsigned char>::iterator>(
        std::vector<unsigned char>::const_iterator &it_in, std::vector<unsigned char>::iterator &it_out);

template
void tidysq::v2::internal::unpack_5<8, std::vector<unsigned char>::const_iterator , std::vector<unsigned char>::iterator>(
        std::vector<unsigned char>::const_iterator &it_in, std::vector<unsigned char>::iterator &it_out);
