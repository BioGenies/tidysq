#pragma once 

#include "tidysq/Sq.h"

namespace tidysq::internal {
        template<typename INTERNAL>
        class MotifFrame {
            std::list<std::string> names_{};
            Sq<INTERNAL> found_;
            std::list<std::string> sought_{};
            std::list<LenSq> start_{};
            std::list<LenSq> end_{};

        public:
            explicit MotifFrame(const Alphabet &alphabet) :
                    found_(Sq<INTERNAL>(alphabet)) {}

            inline void append(const std::string &name,
                               const Sequence<INTERNAL> &found,
                               const std::string &sought,
                               const LenSq start,
                               const LenSq end) {
                names_.push_back(name);
                found_.push_back(found);
                sought_.push_back(sought);
                start_.push_back(start);
                end_.push_back(end);
            }

            inline void merge_with(MotifFrame<INTERNAL> &&other) {
                names_.splice(names_.end(), other.names_);
                sought_.splice(sought_.end(), other.sought_);
                start_.splice(start_.end(), other.start_);
                end_.splice(end_.end(), other.end_);
                for (LenSq i = 0; i < other.found_.size(); i++) {
                    found_.push_back(other.found_[i]);
                }
            }

            [[nodiscard]] const std::list<std::string> &names() const {
                return names_;
            }

            const Sq<INTERNAL> &found() const {
                return found_;
            }

            [[nodiscard]] const std::list<std::string> &sought() const {
                return sought_;
            }

            [[nodiscard]] const std::list<LenSq> &start() const {
                return start_;
            }

            [[nodiscard]] const std::list<LenSq> &end() const {
                return end_;
            }
        };
    }