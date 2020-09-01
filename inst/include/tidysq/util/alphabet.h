#ifndef TIDYSQ_ALPHABETS_H
#define TIDYSQ_ALPHABETS_H

#include <stdexcept>

#include "tidysq/types/Alphabet.h"

namespace tidysq::util {
    template<int DUMMY>
    std::string getDefaultNALetter() {
        return "!";
    }

    template<int DUMMY>
    Alphabet getStandardAlphabet(const SqType &type) {
        std::vector<std::string> letters;
        switch (type) {
            case AMI:
                letters = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                           "S", "T", "U", "V", "W", "X", "Y", "Z", "-", "*"};
                break;
            case AMI_CLN:
                letters = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                           "W", "Y", "-", "*"};
                break;
            case DNA:
                letters = {"A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                break;
            case DNA_CLN:
                letters = {"C", "A", "G", "T", "-"};
                break;
            case RNA:
                letters = {"A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-"};
                break;
            case RNA_CLN:
                letters = {"C", "A", "G", "U", "-"};
                break;
            default:
                throw std::invalid_argument("Provided type does not have a predefined standard alphabet!");
        }
        return Alphabet(letters, getDefaultNALetter<0>());
    }

    template<int DUMMY>
    SqType guessSqType(const Alphabet &alphabet) {
        for (auto &type : {DNA, DNA_CLN, RNA, RNA_CLN, AMI, AMI_CLN}) {
            if (getStandardAlphabet<0>(type) == alphabet) return type;
        }
        return UNT;
    }
}

#endif //TIDYSQ_ALPHABET_H
