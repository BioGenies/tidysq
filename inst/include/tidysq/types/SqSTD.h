#ifndef TIDYSQ_SQSTD_H
#define TIDYSQ_SQSTD_H

#include <utility>
#include <vector>

#include "general.h"
#include "SequenceSTD.h"
#include "Alphabet.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"
#include "tidysq/util/alphabet.h"

namespace tidysq {
    template<>
    class Sq<STD> {
        std::vector<Sequence<STD>> content_;
        Alphabet alphabet_;
        SqType type_;
    public:
        typedef Sequence<STD> SequenceType;

        Sq(const lensq length, const Alphabet &alphabet, const SqType &type) :
                content_(std::vector<Sequence<STD>>(length)),
                alphabet_(alphabet),
                type_(type) {};

        Sq(const Alphabet &alphabet, const SqType &type) :
                Sq(0, alphabet, type) {};

        Sq(const lensq length, const SqType &type) :
                Sq(length, util::getStandardAlphabet(type), type) {};

        explicit Sq(const SqType &type) :
                Sq(util::getStandardAlphabet(type), type) {};

        Sq(const lensq length, const Alphabet &alphabet) :
                Sq(length, alphabet, util::guessSqType(alphabet)) {};

        explicit Sq(const Alphabet& alphabet) :
                Sq(alphabet, util::guessSqType(alphabet)) {};

        SequenceType &operator[] (lensq index) {
            return content_[index];
        }

        const SequenceType &operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] const Alphabet& alphabet() const {
            return alphabet_;
        }

        [[nodiscard]] const SqType& type() const {
            return type_;
        }

        template<InternalType INTERNAL_OUT,
                ProtoType PROTO_OUT>
        SqProto<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<STD>, SqProto<INTERNAL_OUT, PROTO_OUT>>(*this, ops::OperationUnpack<STD, INTERNAL_OUT, PROTO_OUT>());
        }

        void pushBack(const SequenceType &sequence) {
            content_.push_back(sequence);
        }
    };
}

#endif //TIDYSQ_SQSTD_H
