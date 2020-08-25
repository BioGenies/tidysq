#ifndef TIDYSQ_SQSTD_H
#define TIDYSQ_SQSTD_H

#include <utility>
#include <vector>

#include "general.h"
#include "SequenceSTD.h"
#include "Alphabet.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"

namespace tidysq {
    template<>
    class Sq<STD> {
        std::vector<Sequence<STD>> content_;
        Alphabet alphabet_;
    public:
        typedef Sequence<STD> SequenceType;

        explicit Sq(Alphabet alphabet) :
                Sq(0, std::move(alphabet)) {};

        Sq(lensq length, Alphabet alphabet) :
                content_(std::vector<Sequence<STD>>(length)),
                alphabet_(std::move(alphabet)) {};

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

        template<InternalType INTERNAL_OUT,
                ProtoType PROTO_OUT>
        SqProto<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<STD>, SqProto<INTERNAL_OUT, PROTO_OUT>>(*this, ops::OperationUnpack<STD, INTERNAL_OUT, PROTO_OUT>());
        }
    };
}

#endif //TIDYSQ_SQSTD_H
