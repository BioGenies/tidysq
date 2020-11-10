#include <Rcpp.h>
#include <iterator>

#include "tidysq/exports.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_find_invalid_letters(const Rcpp::List& x,
                                    const Rcpp::StringVector& dest_type,
                                    const Rcpp::StringVector& NA_letter) {
    const Sq<RCPP> sq = import_from_R(x, NA_letter);
    const Alphabet& alph = sq.alphabet();
    const std::vector<std::string> dest_alph =
            util::standard_letters_for_sq_type(util::sq_type_for_sq_type_abbr(util::get_scalar_string_value(dest_type)));

    std::vector<LetterValue> invalid_indices;
    for (LetterValue i = 0; i < alph.length(); ++i) {
        if (std::none_of(dest_alph.begin(), dest_alph.end(),
                         [alph, i](const std::string& letter){ return alph[i] == letter; })) {
            invalid_indices.push_back(i);
        }
    }

    Rcpp::List ret = Rcpp::List::create();

    for (LenSq i = 0; i < sq.length(); ++i) {
        const Sequence<RCPP> sequence = sq[i];
        Rcpp::StringVector ret_vector = Rcpp::StringVector();
        std::vector<LetterValue> invalid_found = {};
        for (const LetterValue& index : invalid_indices) {
            if (std::any_of(sequence.cbegin(alph.alphabet_size()), sequence.cend(alph.alphabet_size()),
                    [=](const ElementPacked elem){ return elem == index; })) {
                ret_vector.push_back(alph[index]);
            }
        }
        ret.push_back(ret_vector);
    }

    return ret;
}
