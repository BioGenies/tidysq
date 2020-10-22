#ifndef TIDYSQ_PROTOSEQUENCEINPUTINTERPRETER_H
#define TIDYSQ_PROTOSEQUENCEINPUTINTERPRETER_H

#include "tidysq/types/ProtoSequence.h"

namespace tidysq {
    template<InternalType INTERNAL, ProtoType PROTO, bool SIMPLE>
    class ProtoSequenceInputInterpreter {
        typedef typename ProtoSequence<INTERNAL, PROTO>::ContentType ContentType;
        typedef typename ProtoSequence<INTERNAL, PROTO>::ElementType ElementType;
        typedef typename ContentType::const_iterator ContentConstIteratorType;

              ContentConstIteratorType internal_iterator_;
        const ContentConstIteratorType end_;
        const Alphabet &alphabet_;
        bool reached_end_;

        ProtoSequenceInputInterpreter(const ContentConstIteratorType &iterator_begin,
                                      const ContentConstIteratorType &iterator_end,
                                      const Alphabet &alphabet) :
            internal_iterator_(iterator_begin),
            end_(iterator_end),
            reached_end_(iterator_begin == iterator_end),
            alphabet_(alphabet) {};

    public:
        inline ProtoSequenceInputInterpreter<INTERNAL, PROTO, SIMPLE>& operator++() {
            if (!reached_end_) {
                internal_iterator_++;
                if (internal_iterator_ == end_) reached_end_ = true;
            }
            return *this;
        }

        inline LetterValue operator*() const {
            return reached_end_ ? 0 : alphabet_.match_value((ElementType) *internal_iterator_);
        }

        [[nodiscard]] inline bool reached_end() const {
            return reached_end_;
        }

        friend class ProtoSequence<INTERNAL, PROTO>;
    };
}

#endif //TIDYSQ_PROTOSEQUENCEINPUTINTERPRETER_H
