#pragma once

#include "tidysq/types/internal/LetterTree.h"

namespace tidysq {
    template<InternalType INTERNAL, ProtoType PROTO, bool SIMPLE>
    class ProtoSequenceInputInterpreter {
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType ContentType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceElementType ElementType;
        typedef typename ContentType::const_iterator ContentConstIteratorType;

        ContentConstIteratorType internal_iterator_;
        const ContentConstIteratorType end_;
        const Alphabet &alphabet_;

        bool reached_end_;
        LenSq interpreted_letters_;

        ProtoSequenceInputInterpreter(ContentConstIteratorType &&iterator_begin,
                                      ContentConstIteratorType &&iterator_end,
                                      const Alphabet &alphabet) :
            internal_iterator_(iterator_begin),
            end_(iterator_end),
            reached_end_(internal_iterator_ == end_),
            alphabet_(alphabet),
            interpreted_letters_(0) {};

        [[nodiscard]] inline LetterValue match() const {
            return alphabet_.match_value((ElementType) *internal_iterator_);
        }

    public:
        inline LetterValue get_next() {
            if (reached_end_) {
                return 0;
            } else {
                LetterValue ret = match();
                internal_iterator_++;
                interpreted_letters_++;
                if (internal_iterator_ == end_) reached_end_ = true;
                return ret;
            }
        }

        [[nodiscard]] inline bool reached_end() const {
            return reached_end_;
        }

        [[nodiscard]] inline LenSq interpreted_letters() const {
            return interpreted_letters_;
        }

        friend class ProtoSequence<INTERNAL, PROTO>;
    };

    template<>
    inline LetterValue ProtoSequenceInputInterpreter<RCPP, STRINGS, true>::match() const {
        return alphabet_.match_value((ElementType) internal_iterator_[0]);
    }

    template<>
    inline LetterValue ProtoSequenceInputInterpreter<RCPP, STRINGS, false>::match() const {
        return alphabet_.match_value((ElementType) internal_iterator_[0]);
    }

    template<InternalType INTERNAL>
    class ProtoSequenceInputInterpreter<INTERNAL, STRING, false> {
        typedef typename TypeMapper<INTERNAL, STRING>::ProtoSequenceContentType ContentType;
        typedef ElementStringMultichar ElementType;
        typedef typename ContentType::const_iterator ContentConstIteratorType;

        const Alphabet &alphabet_;
        internal::LetterTree letter_tree_;

        LenSq interpreted_letters_;

        ProtoSequenceInputInterpreter(ContentConstIteratorType &&iterator_begin,
                                      ContentConstIteratorType &&iterator_end,
                                      const Alphabet &alphabet) :
                alphabet_(alphabet),
                letter_tree_(internal::LetterTree(alphabet, std::move(iterator_begin), std::move(iterator_end))),
                interpreted_letters_(0) {
        };

        [[nodiscard]] inline LetterValue match() {
            return letter_tree_.match_next();
        }

    public:
        inline LetterValue get_next() {
            if (reached_end()) {
                return 0;
            } else {
                LetterValue ret = match();
                interpreted_letters_++;
                return ret;
            }
        }

        [[nodiscard]] inline bool reached_end() const {
            return letter_tree_.reached_end();
        }

        [[nodiscard]] inline LenSq interpreted_letters() const {
            return interpreted_letters_;
        }

        friend class ProtoSequence<INTERNAL, STRING>;
    };
}

