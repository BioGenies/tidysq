#pragma once

#include "tidysq/internal/LetterTree.h"

namespace tidysq {
    template<typename INTERNAL, typename PROTO, bool SIMPLE>
    class ProtoSequenceInputInterpreter {
        typedef typename util::TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentStorageType    ContentStorageType;
        typedef typename PROTO::ProtoSequenceElementType                                    ElementType;
        typedef typename ContentStorageType::const_iterator                                 ContentConstIteratorType;

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
            alphabet_(alphabet),
            reached_end_(internal_iterator_ == end_),
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

        inline ElementType get_or_extract_next_element() {
            return get_next_element();
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
    inline LetterValue ProtoSequenceInputInterpreter<RCPP_IT, STRINGS_PT, true>::match_value() const {
        return alphabet_.match_value((ElementType) internal_iterator_[0]);
    }

    template<>
    inline LetterValue ProtoSequenceInputInterpreter<RCPP_IT, STRINGS_PT, false>::match_value() const {
        return alphabet_.match_value((ElementType) internal_iterator_[0]);
    }

    template<typename INTERNAL>
    class ProtoSequenceInputInterpreter<INTERNAL, STRING_PT, false> {
        typedef typename util::TypeBinder<INTERNAL, STRING_PT>::ProtoSequenceContentStorageType ContentStorageType;
        typedef ElementStringMultichar ElementType;
        typedef typename ContentStorageType::const_iterator ContentConstIteratorType;

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

        [[nodiscard]] inline ElementType match_or_extract_element() {
            return letter_tree_.match_or_extract_element();
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

        inline ElementType get_or_extract_next_element() {
            if (reached_end()) {
                return "";
            } else {
                ElementType ret = match_or_extract_element();
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

        friend class ProtoSequence<INTERNAL, STRING_PT>;
    };
}

