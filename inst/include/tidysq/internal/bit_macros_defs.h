// COMMON MACROS:

#define FETCH_OR_ZERO(reg_num, len_condition) \
    if constexpr(LEN > len_condition) { \
        FETCH(reg_num) \
    } else v##reg_num = 0;
#define FETCH(reg_num) \
    v##reg_num = *it_in;\
    ++it_in;
#define FETCH_OR_RETURN(reg_num, len_condition) \
    if constexpr(LEN > len_condition) { \
        FETCH(reg_num) \
    } else return;

// MACROS FOR PACKING:

#define ALIGN_2(reg_num_a, shift_a, reg_num_b, shift_b) \
    *it_out = (v##reg_num_a >> shift_a##u) | (v##reg_num_b << shift_b##u); \
    ++it_out;
#define ALIGN_3(reg_num_a, shift_a, reg_num_b, shift_b, reg_num_c, shift_c) \
    *it_out = (v##reg_num_a >> shift_a##u) | (v##reg_num_b << shift_b##u) | (v##reg_num_c << shift_c##u); \
    ++it_out;
#define ALIGN_4(reg_num_a, shift_a, reg_num_b, shift_b, reg_num_c, shift_c, reg_num_d, shift_d) \
    *it_out = (v##reg_num_a >> shift_a##u) | (v##reg_num_b << shift_b##u) | (v##reg_num_c << shift_c##u) | (v##reg_num_d << shift_d##u); \
    ++it_out;

// MACROS FOR UNPACKING:

#define ALIGN_AND_TRIM_1(reg_num, shift) \
    *it_out = (v##reg_num >> shift##u) & full_byte_trim; \
    ++it_out;
#define ALIGN_AND_TRIM_2(reg_num_a, shift_a, trim_a, reg_num_b, shift_b) \
    *it_out = ((v##reg_num_a >> shift_a##u) & trim_a##u) | ((v##reg_num_b << shift_b##u) & full_byte_trim); \
    ++it_out;
#define ALIGN_AND_TRIM_1_OR_RETURN(reg_num, shift, len_condition) \
    if constexpr(LEN > len_condition) { \
        ALIGN_AND_TRIM_1(reg_num, shift) \
    } else return;
