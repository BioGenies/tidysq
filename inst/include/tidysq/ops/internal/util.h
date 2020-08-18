#ifndef TIDYSQ_UTIL_H
#define TIDYSQ_UTIL_H

#include <cmath>
#include "../../types/SequenceSTD.h"
#include "../../types/SequenceRCPP.h"
#include "../../types/AlphabetRCPP.h"
#include "../../types/AlphabetSTD.h"


namespace tidysq::internal {
    template<InternalType INTERNAL,
            ProtoType PROTO>
    inline lensq getPackedLength(const SequenceProto<INTERNAL, PROTO> &unpacked, const Alphabet<INTERNAL> &alphabet) {
        return (alphabet.alphabetSize() * unpacked.size() + 7) / 8;
    }

    template<InternalType INTERNAL>
    inline lensq getOriginalLength(const Sequence<INTERNAL> &packed);

    template<>
    inline lensq getOriginalLength<STD>(const Sequence<STD> &packed) {
        return packed.originalLength();
    }

    template<>
    inline lensq getOriginalLength<RCPP>(const Sequence<RCPP> &packed) {
        return packed.attr("original_length");
    }

    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserveSpaceForPacked(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                                        const Alphabet<INTERNAL_IN> &alphabet) {
        return Sequence<INTERNAL_OUT>(getPackedLength(unpacked, alphabet), unpacked.size());
    }

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT,
            ProtoType PROTO_OUT>
    inline SequenceProto<INTERNAL_OUT, PROTO_OUT> reserveSpaceForUnpacked(const Sequence<INTERNAL_IN> &packed) {
        return SequenceProto<INTERNAL_OUT, PROTO_OUT>(getOriginalLength(packed));
    }

    template<InternalType INTERNAL_IN>
    struct ValueToLetterMatcher;

    template<>
    struct ValueToLetterMatcher<STD> {
        inline static letvalue match(const std::string &letter, const Alphabet<STD> &alphabet) {
            for (letvalue i = 0; i < alphabet.size(); i++) {
                if (letter == alphabet[i]) {
                    return i + 1;
                }
            }
            return alphabet.NAValue();
        }

        inline static letvalue matchStandard(const char &letter, const Alphabet<STD> &alphabet) {
            for (letvalue i = 0; i < alphabet.size(); i++) {
                if (letter == alphabet[i][0]) {
                    return i + 1;
                }
            }
            return alphabet.NAValue();
        }
    };

    template<>
    struct ValueToLetterMatcher<RCPP> {
        inline static letvalue match(const Rcpp::StringVector::const_Proxy &letter, const Alphabet<RCPP> &alphabet) {
            if (Rcpp::StringVector::is_na(letter)) {
                return alphabet.NAValue();
            }
            for (letvalue i = 0; i < alphabet.size(); i++) {
                if (alphabet[i] == letter) {
                    return i + 1;
                }
            }
            return alphabet.NAValue();
        }
    };

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    struct LetterToValueMatcher;

    template<>
    struct LetterToValueMatcher<RCPP, STD> {
        inline static std::string match(const unsigned char &value, const Alphabet<RCPP> &alphabet) {
            return value == alphabet.NAValue() ? alphabet.NALetter() : alphabet[value - 1];
        }
    };

    template<>
    struct LetterToValueMatcher<STD, STD> {
        inline static std::string match(const unsigned char &value, const Alphabet<STD> &alphabet) {
            return value == alphabet.NAValue() ? alphabet.NALetter() : alphabet[value - 1];
        }

        inline static char matchStandard(const unsigned char &value, const Alphabet<STD> &alphabet) {
            return value == alphabet.NAValue() ? alphabet.NALetter()[0] : alphabet[value - 1][0];
        }
    };

    template<InternalType INTERNAL_IN>
    struct LetterToValueMatcher<INTERNAL_IN, RCPP> {
        inline static Rcpp::String match(const unsigned char &value, const Alphabet<INTERNAL_IN> &alphabet) {
            return value == alphabet.NAValue() ? Rcpp::String(NA_STRING) : alphabet[value - 1];
        }
    };
}

#endif //TIDYSQ_UTIL_H
