#ifndef TIDYSQ_ALPHABETSTD_H
#define TIDYSQ_ALPHABETSTD_H

#include <utility>
#include <vector>
#include <string>
#include <cmath>

#include "general.h"

namespace tidysq {
    template<>
    class Alphabet<STD> : public std::vector<std::string> {
    private:
        typedef std::vector<std::string> BaseType;
        letvalue NAValue_;
        sizealph alphabetSize_;
        std::string NALetter_;
    public:
        Alphabet(const std::vector<std::string>& letters, std::string NALetter) :
                BaseType(letters),
                NALetter_(std::move(NALetter)) {
            alphabetSize_ = ceil(log2((double) letters.size() + 2));
            NAValue_ = pow(2, alphabetSize_) - 1;
        }


        [[nodiscard]] const letvalue &NAValue() const {
            return NAValue_;
        };

        letvalue &NAValue() {
            return NAValue_;
        }

        [[nodiscard]] std::string NALetter() const {
            return NALetter_;
        }

        [[nodiscard]] const sizealph &alphabetSize() const {
            return alphabetSize_;
        }

        [[nodiscard]] bool isSimple() const {
            for (int i = 0; i < size(); i++) {
                if (at(i).size() > 1) {
                    return false;
                }
            }
            return true;
        }
    };
}

#endif //TIDYSQ_ALPHABETSTD_H
