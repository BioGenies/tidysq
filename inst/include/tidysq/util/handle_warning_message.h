#pragma once

#include <string>
#include <Rcpp.h>

#include "tidysq/constants/WarningLevel.h"

namespace tidysq::util {
    inline void handle_warning_message(const std::string &message,
                                       const constants::WarningLevel &warning_level) {
        switch (warning_level) {
            case constants::LVL_SILENT:
                break;
            case constants::LVL_MESSAGE:
                Rcpp::Rcout << message << std::endl;
                break;
            case constants::LVL_WARNING:
                Rcpp::warning(message);
                break;
            case constants::LVL_ERROR:
                Rcpp::stop(message);
            default:
                throw std::invalid_argument("cannot handle warning - invalid warning level!");
        }
    }

    inline void handle_warning_message(const std::string &message,
                                       const std::string &warning_level) {
        return handle_warning_message(message, constants::WARNING_LEVEL_NAMES.at(warning_level));
    }
}