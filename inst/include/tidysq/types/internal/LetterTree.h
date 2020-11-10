#pragma once

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"

namespace tidysq::internal {
    class LetterNode {
        LetterValue value_;
        bool def_init_;
        std::map<const char, LetterNode> further_nodes_;

    public:
        LetterNode() :
                value_(0),
                def_init_(true),
                further_nodes_({}) {};

        explicit LetterNode(const LetterValue NA_value) :
                value_(NA_value),
                def_init_(false),
                further_nodes_({}) {};

        inline LetterNode &match_or_insert(const char &character, const LetterValue &NA_value) {
            LetterNode &ret = further_nodes_[character];
            if (ret.def_init_) ret.value_ = NA_value;
            return ret;
        }

        [[nodiscard]] inline const LetterNode &match(const char &character) const {
            return further_nodes_.at(character);
        }

        [[nodiscard]] const LetterValue &value() const {
            return value_;
        }

        LetterValue &value() {
            return value_;
        }
    };

    class LetterTree {
        typedef std::string                 ContentType;
        typedef ElementStringMultichar      ElementType;
        typedef ContentType::const_iterator ContentIteratorType;

        const Alphabet &alphabet_;
        LetterNode root_;

        ContentIteratorType iterator_;
        const ContentIteratorType end_;

        void put_letter(const Letter &letter, const LetterValue value) {
            auto letter_it = letter.begin();

            LetterNode *node = &root_;

            while(letter_it != letter.end()) {
                node = &(node->match_or_insert(*letter_it, alphabet_.NA_value()));
                letter_it++;
            }

            node->value() = value;
        }

    public:
        explicit LetterTree(const Alphabet &alphabet,
                            ContentIteratorType &&iterator,
                            ContentIteratorType &&end) :
                alphabet_(alphabet),
                root_(LetterNode(alphabet.NA_value())),
                iterator_(iterator),
                end_(end) {
            for (LetterValue i = 0; i < alphabet.length(); i++) {
                put_letter(alphabet[i], i);
            }
            put_letter(alphabet.NA_letter(), alphabet.NA_value());
        };

        LetterValue match_value() {
            const LetterNode* node = &root_;

            while (iterator_ != end_) {
                try {
                    node = &(node->match(*iterator_));
                    iterator_++;
                } catch (const std::out_of_range &e) {
                    return node->value();
                }
            }

            return node->value();
        }

        ElementType match_element() {
            const LetterNode* node = &root_;

            Letter next_element = "";

            while (iterator_ != end_) {
                try {
                    next_element += *iterator_;
                    node = &(node->match(*iterator_));
                    iterator_++;
                } catch (const std::out_of_range &e) {
                    return next_element;
                }
            }

            return next_element;
        }

        [[nodiscard]] inline bool reached_end() const {
            return iterator_ == end_;
        }
    };
}