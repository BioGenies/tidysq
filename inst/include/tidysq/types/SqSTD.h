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

        Sq(lensq length, Alphabet alphabet, const SqType &type) :
                content_(std::vector<Sequence<STD>>(length)),
                alphabet_(std::move(alphabet)),
                type_(type) {};

        Sq(Alphabet alphabet, const SqType &type) :
                Sq(0, std::move(alphabet), type) {};

        Sq(lensq length, const SqType &type) :
                Sq(length, util::getStandardAlphabet<0>(type), type) {};

        explicit Sq(const SqType &type) :
                Sq(util::getStandardAlphabet<0>(type), type) {};

        Sq(lensq length, Alphabet alphabet) :
                Sq(length, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        explicit Sq(Alphabet alphabet) :
                Sq(std::move(alphabet), util::guessSqType<0>(alphabet)) {};

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
    };
}

#endif //TIDYSQ_SQSTD_H
