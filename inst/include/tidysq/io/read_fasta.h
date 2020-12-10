#pragma once

#include <fstream>
#include <stdexcept>

#include "tidysq/constants/io.h"
#include "tidysq/constants/standard_letters.h"
#include "tidysq/Sq.h"

#include "tidysq/ops/pack.h"

namespace tidysq {
    namespace internal {
        template<typename INTERNAL>
        using NamedSqibble = std::tuple<Sq<INTERNAL>, std::vector<std::string>>;

        class FastaSampler {
            enum State {
                AWAITING,
                READING_SEQUENCE,
                READING_NAME,
                VIBING
            };

            class bad_state_exception : public std::exception {
                const char* message_;

            public:
                explicit bad_state_exception(const char* message) :
                        std::exception(),
                        message_(message) {};

                [[nodiscard]] const char * what() const noexcept override {
                    return message_;
                }
            };

            std::ifstream stream_;

            const LenSq sample_size_;
            const Alphabet mock_alphabet_;

            State current_state_;
            char *input_buffer_;

            std::string sequence_buffer_;
            std::set<Letter> letters_;
            unsigned int read_characters_;

            void parse_sequence_buffer() {
                if (sequence_buffer_.empty()) return;

                ProtoSequence<STD_IT, STRING_PT> proto_sequence(sequence_buffer_);

                if (mock_alphabet_.is_simple()) {
                    auto interpreter = proto_sequence.content_interpreter<true>(mock_alphabet_);

                    while (!interpreter.reached_end()) {
                        letters_.insert(Letter{interpreter.get_next_element()});
                    }
                } else {
                    auto interpreter = proto_sequence.content_interpreter<false>(mock_alphabet_);

                    while (!interpreter.reached_end()) {
                        letters_.insert(interpreter.get_next_element());
                    }
                }

                sequence_buffer_.clear();
            }

            inline void handle_letter(const char letter) {
                switch (letter) {
                    case constants::NEW_LINE_CHAR:
                        if (current_state_ == READING_NAME) {
                            current_state_ = READING_SEQUENCE;
                        }
                        break;
                    case constants::NEW_SEQUENCE_CHAR:
                        parse_sequence_buffer();
                        current_state_ = READING_NAME;
                        break;
                    default:
                        if (current_state_ == READING_SEQUENCE) sequence_buffer_ += letter;
                }
            }

            void parse_read_buffer() {
                char const *letter = input_buffer_;
                const auto read_letters = stream_.gcount();

                //for each letter in the buffer - handle it!
                while (letter != input_buffer_ + read_letters) {
                    handle_letter(*letter++);
                }

                // if less letters than buffer may contain was read - end was reached and
                // buffer should be parsed
                if (read_letters < constants::BUFF_SIZE) {
                    parse_sequence_buffer();
                }
            }

        public:
            FastaSampler(const std::string &file_name,
                         const LenSq sample_size,
                         const Letter &NA_letter,
                         const bool ignore_case) :
                    stream_(std::ifstream(file_name)),
                    sample_size_(sample_size),
                    mock_alphabet_(Alphabet({}, tidysq::SqType::UNT, NA_letter, ignore_case)),
                    current_state_(AWAITING),
                    input_buffer_(new char[constants::BUFF_SIZE]),
                    sequence_buffer_({}),
                    letters_({}),
                    read_characters_(0) {
                if (!stream_.is_open()) throw;
            };

            ~FastaSampler() {
                delete[] input_buffer_;
                stream_.close();
            }

            void sample() {
                // while stream is open and no error happened
                while (stream_.good() && read_characters_ < sample_size_) {
                    // read fragment of a file into a buffer and parse it
                    stream_.read(input_buffer_, constants::BUFF_SIZE);
                    parse_read_buffer();
                    read_characters_ += stream_.gcount();
                }

                letters_.erase(mock_alphabet_.NA_letter());
                if (mock_alphabet_.ignores_case()) {
                    std::set<Letter> letters_to_erase = {};
                    for (const auto &letter : letters_) {
                        if (isalpha(letter[0]) && !isupper(letter[0]))
                            letters_to_erase.insert(letter);
                    }
                    for (const auto &letter : letters_to_erase) {
                        letters_.erase(letter);
                    }
                }

                current_state_ = VIBING;
            }

            inline Alphabet alphabet() {
                if (current_state_ != VIBING)
                    throw bad_state_exception("Sampler has not read file contents yet.");
                return Alphabet(tidysq::util::convert_set_to_vector(letters_));
            }

        };

        template<typename INTERNAL>
        class FastaReader {
            enum State {
                AWAITING,
                READING_SEQUENCE,
                READING_NAME,
                VIBING
            };

            class bad_state_exception : public std::exception {
                const char* message_;

            public:
                explicit bad_state_exception(const char* message) :
                        std::exception(),
                        message_(message) {};

                [[nodiscard]] const char * what() const noexcept override {
                    return message_;
                }
            };

            std::ifstream stream_;
            const Alphabet &alphabet_;

            State current_state_;
            char *input_buffer_;
            ProtoSequence<STD_IT, STRING_PT> proto_sequence_buffer;
            std::string name_buffer_;

            Sq<INTERNAL> sq_;
            std::vector<std::string> names_;

            void parse_sequence_buffer() {
                if (proto_sequence_buffer.size() == 0) return;
                sq_.push_back(pack<STD_IT, STRING_PT, INTERNAL>(proto_sequence_buffer, alphabet_));
                proto_sequence_buffer = ProtoSequence<STD_IT, STRING_PT>(0);
            }

            inline void parse_name_buffer() {
                names_.push_back(name_buffer_);
                name_buffer_.clear();
            }

            inline void handle_letter(const char letter) {
                switch (letter) {
                    case constants::NEW_LINE_CHAR:
                        if (current_state_ == READING_NAME) {
                            parse_name_buffer();
                            current_state_ = READING_SEQUENCE;
                        }
                        break;
                    case constants::NEW_SEQUENCE_CHAR:
                        parse_sequence_buffer();
                        current_state_ = READING_NAME;
                        break;
                    default:
                        if (current_state_ == READING_NAME) name_buffer_ += letter;
                        else proto_sequence_buffer += letter;
                }
            }

            void parse_read_buffer() {
                char const *letter = input_buffer_;
                const auto read_letters = stream_.gcount();

                //for each letter in the buffer - handle it!
                while (letter != input_buffer_ + read_letters) {
                    handle_letter(*letter++);
                }

                // if less letters than buffer may contain was read - end was reached and
                // buffer should be parsed
                if (read_letters < constants::BUFF_SIZE) {
                    parse_sequence_buffer();
                }
            }

        public:
            FastaReader(const std::string &file_name, const Alphabet &alphabet) :
                    stream_(std::ifstream(file_name)),
                    alphabet_(alphabet),
                    current_state_(AWAITING),
                    input_buffer_(new char[constants::BUFF_SIZE]),
                    proto_sequence_buffer({}),
                    name_buffer_({}),
                    sq_(Sq<INTERNAL>(0, alphabet)),
                    names_({}) {
                if (!stream_.is_open()) throw std::out_of_range("Out of range!");
            };

            ~FastaReader() {
                delete[] input_buffer_;
                stream_.close();
            }

            void read() {
                // while stream is open and no error happened
                while (stream_.good()) {
                    // read fragment of a file into a buffer and parse it
                    stream_.read(input_buffer_, constants::BUFF_SIZE);
                    parse_read_buffer();
                }

                current_state_ = VIBING;
            }

            inline NamedSqibble<INTERNAL> sqibble() const {
                if (current_state_ != VIBING)
                    throw bad_state_exception("Reader has not read file contents yet.");
                return std::make_tuple(sq_, names_);
            }
        };
    }

    namespace io {
        template<typename INTERNAL>
        internal::NamedSqibble<INTERNAL> read_fasta(const std::string &file_name,
                                                    const Alphabet &alphabet) {
            internal::FastaReader<INTERNAL> reader(file_name, alphabet);
            reader.read();
            return reader.sqibble();
        }

        inline Alphabet sample_fasta(const std::string &file_name,
                                         const LenSq sample_size = constants::BUFF_SIZE,
                                         const Letter &NA_letter = constants::DEFAULT_NA_LETTER,
                                         const bool ignore_case = constants::DEFAULT_IGNORE_CASE) {
            internal::FastaSampler sampler(file_name, sample_size, NA_letter, ignore_case);
            sampler.sample();
            return sampler.alphabet();
        }
    }

}

