#pragma once

#include "tidysq/Alphabet.h"

namespace tidysq::internal {
    class LetterNode {
        LetterValue value_;
        bool end_node_;
        std::map<const char, LetterNode> children_;

    public:
        LetterNode() :
                value_(0),
                end_node_(false),
                children_({}) {};

        inline LetterNode &operator[](const char &character) {
            LetterNode &ret = children_[character];
            return ret;
        }

        [[nodiscard]] inline const LetterNode &at(const char &character) const {
            return children_.at(character);
        }

        inline bool contains(const char &character) const {
            // Replace with contains() after changing to C++20
            return children_.count(character) > 0;
        }

        [[nodiscard]] const LetterValue &value() const {
            return value_;
        }

        [[nodiscard]] const bool &end_node() const {
            return end_node_;
        }

        void make_end_node(LetterValue value) {
            value_ = value;
            end_node_ = true;
        }
    };

    class LetterTree {
        typedef std::string ContentType;
        typedef ElementStringMultichar ElementType;
        typedef ContentType::const_iterator ContentIteratorType;

        const Alphabet &alphabet_;
        LetterNode root_;
        ContentIteratorType iterator_;
        const ContentIteratorType end_;

        void insert(const Letter &letter, const LetterValue value) {
            LetterNode *node = &root_;

            auto letter_it = letter.begin();
            while (letter_it != letter.end()) {
                node = &((*node)[*letter_it]);
                letter_it++;
            }
            node->make_end_node(value);
        }

        [[nodiscard]] inline LetterValue extract_value(const LetterNode &node) const {
            return node.end_node() ? node.value() : alphabet_.NA_value();
        }

        [[nodiscard]] inline Letter extract_element(const LetterNode &node) const {
            return node.end_node() ? alphabet_[node.value()] : alphabet_.NA_letter();
        }

        LetterNode find_next_node() {
            const LetterNode *current_node = &root_;
            const LetterNode *node = &root_;
            ContentIteratorType element_it = iterator_;

            while (element_it != end_) {
                if (current_node->contains(*element_it)) {
                    current_node = &(current_node->at(*element_it));
                    ++element_it;
                    if (current_node->end_node()) {
                        node = current_node;
                        iterator_ = element_it;
                    }
                } else {
                    break;
                }
            }
            return *node;
        }

    public:
        explicit LetterTree(const Alphabet &alphabet,
                            ContentIteratorType &&iterator,
                            ContentIteratorType &&end) :
                alphabet_(alphabet),
                root_(LetterNode()),
                iterator_(iterator),
                end_(end) {
            for (LetterValue i = 0; i < alphabet.size(); i++) {
                insert(alphabet[i], i);
            }
            insert(alphabet.NA_letter(), alphabet.NA_value());
        };

        LetterValue match_value() {
            if (root_.contains(*iterator_)) {
                return extract_value(find_next_node());
            } else {
                ++iterator_;
                return alphabet_.NA_value();
            }
        }

        ElementType match_element() {
            if (root_.contains(*iterator_)) {
                return extract_element(find_next_node());
            } else {
                ++iterator_;
                return alphabet_.NA_letter();
            }
        }

        ElementType match_or_extract_element() {
            // Used for multichar alphabets, when we want to obtain alphabet
            if (root_.contains(*iterator_)) {
                return match_element();
            } else {
                Letter ret = {*iterator_};
                ++iterator_;
                return ret;
            }
        }

        [[nodiscard]] inline bool reached_end() const {
            return iterator_ == end_;
        }
    };
}