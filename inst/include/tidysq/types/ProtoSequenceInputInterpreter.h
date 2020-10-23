#pragma once

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

        struct LetterTree;

        const Alphabet &alphabet_;
        LetterTree letter_tree_;

        LenSq interpreted_letters_;

        ProtoSequenceInputInterpreter(ContentConstIteratorType &&iterator_begin,
                                      ContentConstIteratorType &&iterator_end,
                                      const Alphabet &alphabet) :
                alphabet_(alphabet),
                letter_tree_(LetterTree(alphabet, std::move(iterator_begin), std::move(iterator_end))),
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

    template<InternalType INTERNAL>
    struct ProtoSequenceInputInterpreter<INTERNAL, STRING, false>::LetterTree {
        struct LetterNode;

        const Alphabet &alphabet_;
        const LetterNode &empty_node_;
        std::vector<LetterNode> initial_nodes_;
        ContentConstIteratorType iterator_;
        const ContentConstIteratorType end_;

        void put_letter(const Letter &letter, const LetterValue value) {
            auto next_character = letter.begin();
            auto letter_end = letter.end();

            auto curr_node_list = initial_nodes_;
            auto next_node = curr_node_list.begin();
            auto nodes_list_end = curr_node_list.end();

            while (next_node != nodes_list_end) {
                if (next_node->character_ == *next_character) {
                    next_character++;
                    curr_node_list = next_node->possible_nodes_;
                    next_node = curr_node_list.begin();
                    nodes_list_end = curr_node_list.end();
                } else {
                    next_node++;
                }
            }

            LetterNode &node = *next_node;

            while (next_character != letter_end) {
                node = LetterNode(*next_character, alphabet_.NA_value());
                curr_node_list.push_back(node);
                curr_node_list = node.possible_nodes_;
                next_character++;
            }
            node.value_ = value;
        }

        explicit LetterTree(const Alphabet &alphabet,
                            ContentConstIteratorType &&iterator,
                            ContentConstIteratorType &&end) :
                alphabet_(alphabet),
                empty_node_(LetterNode('\0', alphabet.NA_value())),
                initial_nodes_({}),
                iterator_(iterator),
                end_(end) {
            for (LetterValue i = 0; i < alphabet.length(); i++) {
                put_letter(alphabet[i], i);
            }
            put_letter(alphabet.NA_letter(), alphabet.NA_value());
        };

        LetterValue match_next() {
            std::vector<LetterNode> curr_node_list = initial_nodes_;

            auto next_node = curr_node_list.begin();
            auto nodes_list_end = curr_node_list.end();

            const LetterNode* node = &empty_node_;

            while (curr_node_list.size() > 0 &&
                   next_node != nodes_list_end &&
                   iterator_ != end_) {
                if (next_node->character_ == *iterator_) {
                    node = &(*next_node);
                    curr_node_list = node->possible_nodes_;
                    next_node = curr_node_list.begin();
                    nodes_list_end = curr_node_list.end();
                    iterator_++;
                } else {
                    next_node++;
                }
            }

            return node->value_;
        }

        [[nodiscard]] inline bool reached_end() const {
            return iterator_ == end_;
        }
    };

    template<InternalType INTERNAL>
    struct ProtoSequenceInputInterpreter<INTERNAL, STRING, false>::LetterTree::LetterNode {
        char character_;
        LetterValue value_;
        std::vector<LetterNode> possible_nodes_;

        explicit LetterNode(char character, const LetterValue NA_value) :
                character_(character),
                value_(NA_value),
                possible_nodes_({}) {};
    };


}

