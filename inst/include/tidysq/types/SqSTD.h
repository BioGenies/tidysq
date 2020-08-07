#ifndef TIDYSQ_SQSTD_H
#define TIDYSQ_SQSTD_H

#include <utility>
#include <vector>

#include "general.h"
#include "SequenceSTD.h"
#include "AlphabetSTD.h"
#include "../ops/OperationUnpack.h"
#include "../sqapply.h"

namespace tidysq {
    template<>
    class Sq<STD> {
        std::vector<Sequence<STD>> content_;
        Alphabet<STD> alphabet_;
    public:
        typedef Sequence<STD> SequenceType;
        typedef Alphabet<STD> AlphabetType;

        Sq(lensq length, AlphabetType alphabet) :
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

        [[nodiscard]] AlphabetType alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT,
                ProtoType PROTO_OUT>
        SqProto<INTERNAL_OUT, PROTO_OUT> unpack() const {
            return sqapply<Sq<STD>, SqProto<INTERNAL_OUT, PROTO_OUT>, AlphabetType>(*this, ops::OperationUnpack<STD, INTERNAL_OUT, PROTO_OUT>());
        }
    };
}

#endif //TIDYSQ_SQSTD_H
