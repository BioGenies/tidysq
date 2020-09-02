#ifndef TIDYSQ_READFASTA_H
#define TIDYSQ_READFASTA_H

#include <fstream>
#include <stdexcept>

#include "tidysq/types/all.h"
#include "tidysq/ops/internal/packSTRING.h"

namespace tidysq {
    namespace constants {
        const short BUFF_SIZE = 4096;
        const char NEW_SEQUENCE_CHAR = '>';
        const char NEW_LINE_CHAR = '\n';
    }

    namespace internal {
        inline void parseNameBuffer(std::string &nameBuffer, std::vector<std::string> &names) {
            names.push_back(nameBuffer);
            nameBuffer.clear();
        }

        template<InternalType INTERNAL>
        void parseSequenceBuffer(std::vector<std::string> &sequenceBuffer, Sq<INTERNAL> &sq) {
            if (sequenceBuffer.empty()) return;
            sq.pushBack(
                    packSTRINGS<STD, INTERNAL>(SequenceProto<STD, STRINGS>(sequenceBuffer), sq.alphabet()));
            sequenceBuffer.clear();
        }

        template<InternalType INTERNAL>
        void parseReadBuffer(const char *readBuffer, const unsigned int redCharacters,
                             std::vector<std::string> &sequenceBuffer, std::string &nameBuffer, bool &readingName,
                             Sq<INTERNAL> &sq, std::vector<std::string> &names) {
            char const *letter = readBuffer;
            while (letter != readBuffer + redCharacters) {
                switch (*letter) {
                    case constants::NEW_LINE_CHAR:
                        if (readingName) {
                            parseNameBuffer(nameBuffer, names);
                            readingName = false;
                        }
                        break;
                    case constants::NEW_SEQUENCE_CHAR:
                        parseSequenceBuffer(sequenceBuffer, sq);
                        readingName = true;
                        break;
                    default:
                        if (readingName) nameBuffer += *letter;
                        else sequenceBuffer.emplace_back(std::string({*letter}));
                }
                letter++;
            }
            if (redCharacters < constants::BUFF_SIZE) {
                parseSequenceBuffer(sequenceBuffer, sq);
            }
        }

        inline void parseReadBuffer(const char *readBuffer, const unsigned int redCharacters,
                                    std::set<char> &alphabetBuffer, bool &readingName) {
            char const *letter = readBuffer;
            while (letter != readBuffer + redCharacters) {
                switch (*letter) {
                    case constants::NEW_LINE_CHAR:
                        if (readingName) {
                            readingName = false;
                        }
                        break;
                    case constants::NEW_SEQUENCE_CHAR:
                        readingName = true;
                        break;
                    default:
                        if (!readingName) alphabetBuffer.emplace(*letter);
                }
            }
        }

        inline Alphabet sampleFasta(const std::string &fileName, const unsigned int sampleSize) {
            std::ifstream stream(fileName);

            if (!stream.is_open()) throw;

            char *readBuffer = new char[constants::BUFF_SIZE];
            bool readingName = false;
            std::set<char> alphabetBuffer{};
            unsigned int leftToSample = sampleSize;
            unsigned int redCharacters = 0;

            while (stream.good() && leftToSample > 0) {
                stream.read(readBuffer, constants::BUFF_SIZE);
                redCharacters = stream.gcount();
                if (redCharacters > leftToSample) redCharacters = leftToSample;

                internal::parseReadBuffer(readBuffer, redCharacters, alphabetBuffer, readingName);

                leftToSample -= redCharacters;
            }

            delete[] readBuffer;
            stream.close();

            return Alphabet(alphabetBuffer, util::getDefaultNALetter());
        }
    }

    template<InternalType INTERNAL>
    std::tuple<Sq<INTERNAL>, std::vector<std::string>> readFasta(const std::string &fileName, const Alphabet &alphabet,
                                                                 const SqType &type) {
        std::ifstream stream(fileName);
        if (!stream.is_open()) throw;

        char *readBuffer = new char[constants::BUFF_SIZE];
        std::vector<std::string> sequenceBuffer{};
        std::string nameBuffer{};
        bool readingName = false;
        Sq<INTERNAL> sq(0, alphabet, type);
        std::vector<std::string> names{};

        while (stream.good()) {
            stream.read(readBuffer, constants::BUFF_SIZE);
            internal::parseReadBuffer(readBuffer, stream.gcount(), sequenceBuffer, nameBuffer, readingName, sq, names);
        }

        delete[] readBuffer;
        stream.close();

        return std::make_tuple(sq, names);
    }

    template<InternalType INTERNAL>
    std::tuple<Sq<INTERNAL>, std::vector<std::string>> readFasta(const std::string &fileName, const Alphabet &alphabet) {
        return readFasta<INTERNAL>(fileName, alphabet, util::guessSqType(alphabet));
    }

    template<InternalType INTERNAL>
    std::tuple<Sq<INTERNAL>, std::vector<std::string>> readFasta(const std::string &fileName, const SqType &type) {
        return readFasta<INTERNAL>(fileName, util::getStandardAlphabet(type), type);
    }

    template<InternalType INTERNAL>
    std::tuple<Sq<INTERNAL>, std::vector<std::string>> readFasta(const std::string &fileName, const unsigned int sampleSize = 4096) {
        Alphabet alphabet = internal::sampleFasta(fileName, sampleSize);
        return readFasta<INTERNAL>(fileName, alphabet);
    }
}

#endif //TIDYSQ_READFASTA_H
