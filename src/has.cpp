#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

std::list<std::list<std::list<std::string>>> convert_motifs(const Rcpp::StringVector& motifs,
                                                            const Alphabet& alph) {
    std::list<std::list<std::list<std::string>>> ret{};
    for (const auto& motif : motifs) {
        auto motif_string = Rcpp::as<std::string>(motif);
        std::list<std::list<std::string>> motif_as_list{};
        for (const char& letter : motif_string) {
            motif_as_list.push_back({std::string(1, letter)});
        }
        ret.push_back(motif_as_list);
    }
    return ret;
}

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List& x,
                            const Rcpp::StringVector& motifs) {
    const Sq<RCPP> sq = importFromR(x, "!");
    const Alphabet& alph = sq.alphabet();
    Rcpp::LogicalVector ret(sq.length());

    // TODO: implement possibility of reading motifs for multiletter alphabets
    if (!alph.is_simple()) throw std::exception();
    const std::list<std::list<std::list<std::string>>> motif_list = convert_motifs(motifs, alph);
    for (LenSq i = 0; i < sq.length(); ++i) {
        const Sequence<RCPP> sequence = sq[i];
        ret[i] = std::all_of(motif_list.begin(), motif_list.end(), [=](const std::list<std::list<std::string>>& motif) {
            SequenceIterator<RCPP> it = sequence.begin(alph);
            bool contains_motif = motif.empty();
            if (sequence.originalLength() >= motif.size()) {
                while (!contains_motif && it < sequence.end(alph) - motif.size()) {
                    auto motif_it = motif.begin();
                    SequenceIterator<RCPP> it2 = it;
                    while (it2 <= sequence.end(alph) && std::any_of(
                            motif_it->begin(), motif_it->end(), [=](const std::string& possible_letter) {
                                return alph[*it2] == possible_letter;
                            })) {
                        ++motif_it;
                        ++it2;
                        if (motif_it == motif.end()) {
                            contains_motif = true;
                            break;
                        }
                    }
                    ++it;
                }
            }
            return contains_motif;
        });
    }
    // TODO: die again
    // Steps to take:
    // 1. handle ^ and $ signs
    // 2. convert motifs to a list of lists of letters
    //      (if you go for a class, place it within tidysq-internal namespace)
    //      (each element of inner lists is a set of acceptable letters at this position - stored as bits)
    //      (store one motif as linked list of arrays)
    // 3. iterate over sequences and motifs, trying to align them
    //      (probably for each vector in a list within the list use any_of(letter_vec...))
    return ret;
}
