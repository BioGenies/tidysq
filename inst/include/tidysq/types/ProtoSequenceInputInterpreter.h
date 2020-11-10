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

        [[nodiscard]] inline ElementType match_element() const {
            return (ElementType) *internal_iterator_;
        }

        [[nodiscard]] inline LetterValue match_value() const {
            return alphabet_.match_value(match_element());
        }

    public:
        inline LetterValue get_next_value() {
            if (reached_end_) {
                return 0;
            } else {
                LetterValue ret = match_value();
                internal_iterator_++;
                interpreted_letters_++;
                if (internal_iterator_ == end_) reached_end_ = true;
                return ret;
            }
        }


        inline ElementType get_next_element() {
            if (reached_end_) {
                return (ElementType) 0;
            } else {
                ElementType ret = match_element();
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
    inline LetterValue ProtoSequenceInputInterpreter<RCPP, STRINGS, true>::match_value() const {
        return alphabet_.match_value((ElementType) internal_iterator_[0]);
    }

    template<>
    inline LetterValue ProtoSequenceInputInterpreter<RCPP, STRINGS, false>::match_value() const {
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

        [[nodiscard]] inline ElementType match_element() {
            return letter_tree_.match_element();
        }

        [[nodiscard]] inline LetterValue match_value() {
            return letter_tree_.match_value();
        }

    public:
        inline LetterValue get_next_value() {
            if (reached_end()) {
                return 0;
            } else {
                LetterValue ret = match_value();
                interpreted_letters_++;
                return ret;
            }
        }

        inline ElementType get_next_element() {
            if (reached_end()) {
                return "";
            } else {
                ElementType ret = match_element();
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

