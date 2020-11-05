#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_obtain_alphabet(const Rcpp::StringVector &x,
                                       const Rcpp::IntegerVector &sample_size,
                                       const Rcpp::StringVector &NA_letter) {
    LenSq already_read = 0;
    std::set<Letter> letters = {};
    const Letter NA_letter_s = util::get_scalar_string_value(NA_letter);
    const unsigned int sample_size_s = sample_size[0];

    auto iter = x.begin();

    // TODO: generify it!
    if (NA_letter_s.length() == 1) {
        while(already_read < sample_size_s && iter != x.end()) {
            ProtoSequence<RCPP, STRING> sequence((char *) *iter);
            auto interpreter = sequence.content_interpreter<true>(Alphabet({}, NA_letter_s));
            while (already_read < sample_size_s && !interpreter.reached_end()) {
                letters.insert({interpreter.get_next_element()});
                already_read += 1;
            }
            iter++;
        }
    } else {
        while(already_read < sample_size[0] && iter != x.end()) {
            ProtoSequence<RCPP, STRING> sequence((char *) *iter);
            auto interpreter = sequence.content_interpreter<false>(Alphabet({}, NA_letter_s));
            while (already_read < sample_size_s && !interpreter.reached_end()) {
                letters.insert(interpreter.get_next_element());
                already_read += 1;
            }
            iter++;
        }
    }

    letters.erase(NA_letter_s);

    return export_to_R(Alphabet(util::convert_set_to_vector(letters), UNT, NA_letter_s));
}

