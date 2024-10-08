# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CPP_apply_R_function <- function(x, fun, single_string, NA_letter) {
    .Call(`_tidysq_CPP_apply_R_function`, x, fun, single_string, NA_letter)
}

CPP_bite <- function(x, indices, NA_letter, on_warning) {
    .Call(`_tidysq_CPP_bite`, x, indices, NA_letter, on_warning)
}

CPP_collapse <- function(x, NA_letter) {
    .Call(`_tidysq_CPP_collapse`, x, NA_letter)
}

CPP_complement <- function(x, NA_letter) {
    .Call(`_tidysq_CPP_complement`, x, NA_letter)
}

CPP_find_invalid_letters <- function(x, dest_type, NA_letter) {
    .Call(`_tidysq_CPP_find_invalid_letters`, x, dest_type, NA_letter)
}

CPP_find_motifs <- function(x, names, motifs, NA_letter) {
    .Call(`_tidysq_CPP_find_motifs`, x, names, motifs, NA_letter)
}

CPP_get_standard_alphabet <- function(dest_type) {
    .Call(`_tidysq_CPP_get_standard_alphabet`, dest_type)
}

CPP_guess_standard_alph <- function(alph, NA_letter) {
    .Call(`_tidysq_CPP_guess_standard_alph`, alph, NA_letter)
}

CPP_has <- function(x, motifs, NA_letter) {
    .Call(`_tidysq_CPP_has`, x, motifs, NA_letter)
}

CPP_obtain_alphabet <- function(x, sample_size, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_obtain_alphabet`, x, sample_size, NA_letter, ignore_case)
}

CPP_pack_RAWS <- function(proto, alphabet, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_pack_RAWS`, proto, alphabet, NA_letter, ignore_case)
}

CPP_pack_INTS <- function(proto, alphabet, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_pack_INTS`, proto, alphabet, NA_letter, ignore_case)
}

CPP_pack_STRINGS <- function(proto, alphabet, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_pack_STRINGS`, proto, alphabet, NA_letter, ignore_case)
}

CPP_pack_STRING <- function(proto, alphabet, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_pack_STRING`, proto, alphabet, NA_letter, ignore_case)
}

CPP_paste <- function(list_of_x, NA_letter) {
    .Call(`_tidysq_CPP_paste`, list_of_x, NA_letter)
}

CPP_random_sq <- function(n, len, alphabet, use_gap) {
    .Call(`_tidysq_CPP_random_sq`, n, len, alphabet, use_gap)
}

CPP_read_fasta <- function(file_name, alphabet, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_read_fasta`, file_name, alphabet, NA_letter, ignore_case)
}

CPP_sample_fasta <- function(file_name, sample_size, NA_letter, ignore_case) {
    .Call(`_tidysq_CPP_sample_fasta`, file_name, sample_size, NA_letter, ignore_case)
}

CPP_remove_NA <- function(x, by_letter, NA_letter) {
    .Call(`_tidysq_CPP_remove_NA`, x, by_letter, NA_letter)
}

CPP_remove_ambiguous <- function(x, by_letter, NA_letter) {
    .Call(`_tidysq_CPP_remove_ambiguous`, x, by_letter, NA_letter)
}

CPP_reverse <- function(x, NA_letter) {
    .Call(`_tidysq_CPP_reverse`, x, NA_letter)
}

CPP_substitute_letters <- function(x, encoding, NA_letter) {
    .Call(`_tidysq_CPP_substitute_letters`, x, encoding, NA_letter)
}

CPP_translate <- function(x, table, NA_letter) {
    .Call(`_tidysq_CPP_translate`, x, table, NA_letter)
}

CPP_typify <- function(x, dest_type, NA_letter) {
    .Call(`_tidysq_CPP_typify`, x, dest_type, NA_letter)
}

CPP_unpack_RAWS <- function(sq, NA_letter) {
    .Call(`_tidysq_CPP_unpack_RAWS`, sq, NA_letter)
}

CPP_unpack_INTS <- function(sq, NA_letter) {
    .Call(`_tidysq_CPP_unpack_INTS`, sq, NA_letter)
}

CPP_unpack_STRINGS <- function(sq, NA_letter) {
    .Call(`_tidysq_CPP_unpack_STRINGS`, sq, NA_letter)
}

CPP_unpack_STRING <- function(sq, NA_letter) {
    .Call(`_tidysq_CPP_unpack_STRING`, sq, NA_letter)
}

CPP_write_fasta <- function(x, names, file, width, NA_value) {
    invisible(.Call(`_tidysq_CPP_write_fasta`, x, names, file, width, NA_value))
}

