#ifndef TIDYSQ_SQPROTOSTD_H
#define TIDYSQ_SQPROTOSTD_H

#include <utility>

#include "SequenceProtoSTD.h"
#include "Alphabet.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"
#include "tidysq/util/alphabet.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<STD, PROTO> {
        std::vector<SequenceProto<STD, PROTO>> content_;
        Alphabet alphabet_;
        SqType type_;
    public:
        typedef SequenceProto<STD, PROTO> SequenceType;

        SqProto(const std::vector<SequenceType> &content, Alphabet alphabet, const SqType &type) :
                content_(content),
                alphabet_(std::move(alphabet)),
                type_(type) {};

        SqProto(const lensq length, Alphabet alphabet, const SqType &type) :
                SqProto(std::vector<SequenceType>(length), std::move(alphabet), type) {};

        SqProto(const std::vector<SequenceType> &content, Alphabet alphabet) :
                SqProto(content, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        SqProto(const lensq length, Alphabet alphabet) :
                SqProto(length, std::move(alphabet), util::guessSqType<0>(alphabet)) {};

        SqProto(const std::vector<SequenceType> &content, const SqType &type) :
                SqProto(length, util::getStandardAlphabet<0>(type), type) {};

        SqProto(const lensq length, const SqType &type) :
                SqProto(length, util::getStandardAlphabet<0>(type), type) {};

        SequenceType &operator[] (lensq index) {
            return content_[index];
        }

        const SequenceType &operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] const Alphabet &alphabet() const {
            return alphabet_;
        }

        [[nodiscard]] const SqType &type() const {
            return type_;
        };

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<STD, PROTO>, Sq<INTERNAL_OUT>>(*this, ops::OperationPack<STD, PROTO, INTERNAL_OUT>());
        }
    };
}

#endif //TIDYSQ_SQPROTOSTD_H
