// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/tidysq.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// C_pack_raws
Rcpp::RawVector C_pack_raws(Rcpp::RawVector unpacked, const unsigned short alph_size);
static SEXP _tidysq_C_pack_raws_try(SEXP unpackedSEXP, SEXP alph_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type unpacked(unpackedSEXP);
    Rcpp::traits::input_parameter< const unsigned short >::type alph_size(alph_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_pack_raws(unpacked, alph_size));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_pack_raws(SEXP unpackedSEXP, SEXP alph_sizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_pack_raws_try(unpackedSEXP, alph_sizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_pack_ints
Rcpp::RawVector C_pack_ints(Rcpp::IntegerVector unpacked, const unsigned short alph_size);
static SEXP _tidysq_C_pack_ints_try(SEXP unpackedSEXP, SEXP alph_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type unpacked(unpackedSEXP);
    Rcpp::traits::input_parameter< const unsigned short >::type alph_size(alph_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_pack_ints(unpacked, alph_size));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_pack_ints(SEXP unpackedSEXP, SEXP alph_sizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_pack_ints_try(unpackedSEXP, alph_sizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_pack_chars
Rcpp::RawVector C_pack_chars(Rcpp::CharacterVector unpacked, Rcpp::CharacterVector alph);
static SEXP _tidysq_C_pack_chars_try(SEXP unpackedSEXP, SEXP alphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type unpacked(unpackedSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    rcpp_result_gen = Rcpp::wrap(C_pack_chars(unpacked, alph));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_pack_chars(SEXP unpackedSEXP, SEXP alphSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_pack_chars_try(unpackedSEXP, alphSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_pack_string
Rcpp::RawVector C_pack_string(Rcpp::RawVector unpacked, Rcpp::CharacterVector alph);
static SEXP _tidysq_C_pack_string_try(SEXP unpackedSEXP, SEXP alphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type unpacked(unpackedSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    rcpp_result_gen = Rcpp::wrap(C_pack_string(unpacked, alph));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_pack_string(SEXP unpackedSEXP, SEXP alphSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_pack_string_try(unpackedSEXP, alphSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// nc_pack_cdna
Rcpp::RawVector nc_pack_cdna(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_cdna(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_cdna(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// nc_pack_crna
Rcpp::RawVector nc_pack_crna(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_crna(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_crna(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// nc_pack_dna
Rcpp::RawVector nc_pack_dna(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_dna(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_dna(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// nc_pack_rna
Rcpp::RawVector nc_pack_rna(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_rna(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_rna(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// nc_pack_cami
Rcpp::RawVector nc_pack_cami(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_cami(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_cami(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// nc_pack_ami
Rcpp::RawVector nc_pack_ami(Rcpp::RawVector UNPACKED);
RcppExport SEXP _tidysq_nc_pack_ami(SEXP UNPACKEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type UNPACKED(UNPACKEDSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_pack_ami(UNPACKED));
    return rcpp_result_gen;
END_RCPP
}
// C_unpack_sq_parallel
Rcpp::List C_unpack_sq_parallel(Rcpp::List sq);
RcppExport SEXP _tidysq_C_unpack_sq_parallel(SEXP sqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sq(sqSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unpack_sq_parallel(sq));
    return rcpp_result_gen;
END_RCPP
}
// nc_read_fasta_file
Rcpp::List nc_read_fasta_file(std::string file, std::string type, bool is_clean);
RcppExport SEXP _tidysq_nc_read_fasta_file(SEXP fileSEXP, SEXP typeSEXP, SEXP is_cleanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type is_clean(is_cleanSEXP);
    rcpp_result_gen = Rcpp::wrap(nc_read_fasta_file(file, type, is_clean));
    return rcpp_result_gen;
END_RCPP
}
// read_fasta_file
Rcpp::List read_fasta_file(std::string file, Rcpp::CharacterVector alph);
RcppExport SEXP _tidysq_read_fasta_file(SEXP fileSEXP, SEXP alphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    rcpp_result_gen = Rcpp::wrap(read_fasta_file(file, alph));
    return rcpp_result_gen;
END_RCPP
}
// find_alph
std::list<char> find_alph(std::string file);
RcppExport SEXP _tidysq_find_alph(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(find_alph(file));
    return rcpp_result_gen;
END_RCPP
}
// tmpPack
List tmpPack(List raws, StringVector alphabet);
RcppExport SEXP _tidysq_tmpPack(SEXP rawsSEXP, SEXP alphabetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type raws(rawsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type alphabet(alphabetSEXP);
    rcpp_result_gen = Rcpp::wrap(tmpPack(raws, alphabet));
    return rcpp_result_gen;
END_RCPP
}
// tmpUnpack
List tmpUnpack(List raws);
RcppExport SEXP _tidysq_tmpUnpack(SEXP rawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type raws(rawsSEXP);
    rcpp_result_gen = Rcpp::wrap(tmpUnpack(raws));
    return rcpp_result_gen;
END_RCPP
}
// tmpUnpack2
List tmpUnpack2(List raws);
RcppExport SEXP _tidysq_tmpUnpack2(SEXP rawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type raws(rawsSEXP);
    rcpp_result_gen = Rcpp::wrap(tmpUnpack2(raws));
    return rcpp_result_gen;
END_RCPP
}
// tmpUnpack3
List tmpUnpack3(List raws);
RcppExport SEXP _tidysq_tmpUnpack3(SEXP rawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type raws(rawsSEXP);
    rcpp_result_gen = Rcpp::wrap(tmpUnpack3(raws));
    return rcpp_result_gen;
END_RCPP
}
// C_unpack_raws
Rcpp::RawVector C_unpack_raws(Rcpp::RawVector packed, const unsigned short alph_size);
static SEXP _tidysq_C_unpack_raws_try(SEXP packedSEXP, SEXP alph_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type packed(packedSEXP);
    Rcpp::traits::input_parameter< const unsigned short >::type alph_size(alph_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unpack_raws(packed, alph_size));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_unpack_raws(SEXP packedSEXP, SEXP alph_sizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_unpack_raws_try(packedSEXP, alph_sizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_unpack_ints
Rcpp::IntegerVector C_unpack_ints(Rcpp::RawVector packed, const unsigned short alph_size);
static SEXP _tidysq_C_unpack_ints_try(SEXP packedSEXP, SEXP alph_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type packed(packedSEXP);
    Rcpp::traits::input_parameter< const unsigned short >::type alph_size(alph_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unpack_ints(packed, alph_size));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_unpack_ints(SEXP packedSEXP, SEXP alph_sizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_unpack_ints_try(packedSEXP, alph_sizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_unpack_chars
Rcpp::CharacterVector C_unpack_chars(Rcpp::RawVector packed, Rcpp::CharacterVector alph, Rcpp::CharacterVector na_char);
static SEXP _tidysq_C_unpack_chars_try(SEXP packedSEXP, SEXP alphSEXP, SEXP na_charSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type packed(packedSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type na_char(na_charSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unpack_chars(packed, alph, na_char));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_unpack_chars(SEXP packedSEXP, SEXP alphSEXP, SEXP na_charSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_unpack_chars_try(packedSEXP, alphSEXP, na_charSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_unpack_string
Rcpp::CharacterVector C_unpack_string(Rcpp::RawVector packed, Rcpp::CharacterVector alph, Rcpp::CharacterVector na_char);
static SEXP _tidysq_C_unpack_string_try(SEXP packedSEXP, SEXP alphSEXP, SEXP na_charSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type packed(packedSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type na_char(na_charSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unpack_string(packed, alph, na_char));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_unpack_string(SEXP packedSEXP, SEXP alphSEXP, SEXP na_charSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_unpack_string_try(packedSEXP, alphSEXP, na_charSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_get_alph_size
unsigned short C_get_alph_size(Rcpp::CharacterVector alph);
static SEXP _tidysq_C_get_alph_size_try(SEXP alphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type alph(alphSEXP);
    rcpp_result_gen = Rcpp::wrap(C_get_alph_size(alph));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_get_alph_size(SEXP alphSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_get_alph_size_try(alphSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// C_get_na_val
unsigned short C_get_na_val(const unsigned short alph_size);
static SEXP _tidysq_C_get_na_val_try(SEXP alph_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const unsigned short >::type alph_size(alph_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_get_na_val(alph_size));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidysq_C_get_na_val(SEXP alph_sizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidysq_C_get_na_val_try(alph_sizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _tidysq_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::RawVector(*C_pack_raws)(Rcpp::RawVector,const unsigned short)");
        signatures.insert("Rcpp::RawVector(*C_pack_ints)(Rcpp::IntegerVector,const unsigned short)");
        signatures.insert("Rcpp::RawVector(*C_pack_chars)(Rcpp::CharacterVector,Rcpp::CharacterVector)");
        signatures.insert("Rcpp::RawVector(*C_pack_string)(Rcpp::RawVector,Rcpp::CharacterVector)");
        signatures.insert("Rcpp::RawVector(*C_unpack_raws)(Rcpp::RawVector,const unsigned short)");
        signatures.insert("Rcpp::IntegerVector(*C_unpack_ints)(Rcpp::RawVector,const unsigned short)");
        signatures.insert("Rcpp::CharacterVector(*C_unpack_chars)(Rcpp::RawVector,Rcpp::CharacterVector,Rcpp::CharacterVector)");
        signatures.insert("Rcpp::CharacterVector(*C_unpack_string)(Rcpp::RawVector,Rcpp::CharacterVector,Rcpp::CharacterVector)");
        signatures.insert("unsigned short(*C_get_alph_size)(Rcpp::CharacterVector)");
        signatures.insert("unsigned short(*C_get_na_val)(const unsigned short)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _tidysq_RcppExport_registerCCallable() { 
    R_RegisterCCallable("tidysq", "_tidysq_C_pack_raws", (DL_FUNC)_tidysq_C_pack_raws_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_pack_ints", (DL_FUNC)_tidysq_C_pack_ints_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_pack_chars", (DL_FUNC)_tidysq_C_pack_chars_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_pack_string", (DL_FUNC)_tidysq_C_pack_string_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_unpack_raws", (DL_FUNC)_tidysq_C_unpack_raws_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_unpack_ints", (DL_FUNC)_tidysq_C_unpack_ints_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_unpack_chars", (DL_FUNC)_tidysq_C_unpack_chars_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_unpack_string", (DL_FUNC)_tidysq_C_unpack_string_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_get_alph_size", (DL_FUNC)_tidysq_C_get_alph_size_try);
    R_RegisterCCallable("tidysq", "_tidysq_C_get_na_val", (DL_FUNC)_tidysq_C_get_na_val_try);
    R_RegisterCCallable("tidysq", "_tidysq_RcppExport_validate", (DL_FUNC)_tidysq_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_tidysq_C_pack_raws", (DL_FUNC) &_tidysq_C_pack_raws, 2},
    {"_tidysq_C_pack_ints", (DL_FUNC) &_tidysq_C_pack_ints, 2},
    {"_tidysq_C_pack_chars", (DL_FUNC) &_tidysq_C_pack_chars, 2},
    {"_tidysq_C_pack_string", (DL_FUNC) &_tidysq_C_pack_string, 2},
    {"_tidysq_nc_pack_cdna", (DL_FUNC) &_tidysq_nc_pack_cdna, 1},
    {"_tidysq_nc_pack_crna", (DL_FUNC) &_tidysq_nc_pack_crna, 1},
    {"_tidysq_nc_pack_dna", (DL_FUNC) &_tidysq_nc_pack_dna, 1},
    {"_tidysq_nc_pack_rna", (DL_FUNC) &_tidysq_nc_pack_rna, 1},
    {"_tidysq_nc_pack_cami", (DL_FUNC) &_tidysq_nc_pack_cami, 1},
    {"_tidysq_nc_pack_ami", (DL_FUNC) &_tidysq_nc_pack_ami, 1},
    {"_tidysq_C_unpack_sq_parallel", (DL_FUNC) &_tidysq_C_unpack_sq_parallel, 1},
    {"_tidysq_nc_read_fasta_file", (DL_FUNC) &_tidysq_nc_read_fasta_file, 3},
    {"_tidysq_read_fasta_file", (DL_FUNC) &_tidysq_read_fasta_file, 2},
    {"_tidysq_find_alph", (DL_FUNC) &_tidysq_find_alph, 1},
    {"_tidysq_tmpPack", (DL_FUNC) &_tidysq_tmpPack, 2},
    {"_tidysq_tmpUnpack", (DL_FUNC) &_tidysq_tmpUnpack, 1},
    {"_tidysq_tmpUnpack2", (DL_FUNC) &_tidysq_tmpUnpack2, 1},
    {"_tidysq_tmpUnpack3", (DL_FUNC) &_tidysq_tmpUnpack3, 1},
    {"_tidysq_C_unpack_raws", (DL_FUNC) &_tidysq_C_unpack_raws, 2},
    {"_tidysq_C_unpack_ints", (DL_FUNC) &_tidysq_C_unpack_ints, 2},
    {"_tidysq_C_unpack_chars", (DL_FUNC) &_tidysq_C_unpack_chars, 3},
    {"_tidysq_C_unpack_string", (DL_FUNC) &_tidysq_C_unpack_string, 3},
    {"_tidysq_C_get_alph_size", (DL_FUNC) &_tidysq_C_get_alph_size, 1},
    {"_tidysq_C_get_na_val", (DL_FUNC) &_tidysq_C_get_na_val, 1},
    {"_tidysq_RcppExport_registerCCallable", (DL_FUNC) &_tidysq_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_tidysq(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
