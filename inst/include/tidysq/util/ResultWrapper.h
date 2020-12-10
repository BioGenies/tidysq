#pragma once

#include <optional>
#include <string>

namespace tidysq::util {
    template<typename RESULT_TYPE>
    class ResultWrapper {
    public:
        typedef RESULT_TYPE ResultType;
        typedef std::optional<std::string> MessageType;
    private:
        const ResultType result_;
        const MessageType message_;
    public:
        ResultWrapper(ResultType result,
                      const MessageType &message) :
                result_(std::move(result)),
                message_(message) {
        };

        operator ResultType() const {
            return result();
        }

        [[nodiscard]] ResultType result() const {
            return result_;
        }

        [[nodiscard]] bool has_message() const {
            return message_.has_value();
        }

        [[nodiscard]] std::string message_text() const {
            return *message_;
        }
    };
}