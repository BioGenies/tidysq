#pragma once

#include <fstream>

#include "tidysq/Sq.h"
#include "tidysq/constants/io.h"
#include "tidysq/ops/unpack.h"

namespace tidysq {
    namespace internal {
        template<typename INTERNAL>
        class FastaWriter {
            std::ofstream stream_;
            const unsigned int width_;

            const Sq<INTERNAL> &sq_;
            const std::vector<std::string> &names_;

            void write_name(LenSq i) {
                stream_ << ">" << names_[i] << std::endl;
            }

            void write_sequence_part(const std::string &content,
                                     u_LenSq &written) {
                // if there is more to be written than content size, write only part of it
                if (content.size() - written >= width_) {
                    stream_.write(content.data() + written, width_);
                    written += width_;
                } else { // else - write everything that's left
                    stream_.write(content.data() + written, content.size() - written);
                    written = content.size();
                }
                stream_ << std::endl;
            }

            void write_sequence(LenSq i) {
                ProtoSequence<STD_IT, STRING_PT> unpacked = unpack<INTERNAL, STD_IT, STRING_PT>(sq_[i], sq_.alphabet());
                const std::string &content = unpacked.content();
                u_LenSq written = 0;

                while (written < content.size()) {
                    write_sequence_part(content, written);
                }
            }

        public:
            FastaWriter(const Sq<INTERNAL> &sq,
                        const std::vector<std::string> &names,
                        const std::string &file_name,
                        const unsigned int &width) :
                    stream_(std::ofstream(file_name)),
                    width_(width),
                    sq_(sq),
                    names_(names) {
                if (!stream_.is_open()) throw std::out_of_range("Out of range!");
            };

            ~FastaWriter() {
                stream_.close();
            }

            void write() {
                for (LenSq i = 0; i < sq_.size(); i++) {
                    write_name(i);
                    write_sequence(i);
                }
            }
        };
    }

    namespace io {
        template<typename INTERNAL>
        void write_fasta(const Sq<INTERNAL> &sq,
                         const std::vector<std::string> &names,
                         const std::string &file_name,
                         const unsigned int &width) {
            auto writer = internal::FastaWriter<INTERNAL>(sq, names, file_name, width);
            writer.write();
        }
    }

}