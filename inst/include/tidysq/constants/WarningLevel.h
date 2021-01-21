#pragma once

#include <unordered_map>
#include <string>

namespace tidysq::constants {
    enum WarningLevel {
        LVL_SILENT,
        LVL_MESSAGE,
        LVL_WARNING,
        LVL_ERROR
    };

    const WarningLevel DEFAULT_WARNING_LEVEL = WarningLevel::LVL_WARNING;

    const std::unordered_map<std::string, WarningLevel> WARNING_LEVEL_NAMES {
            {"silent", LVL_SILENT},
            {"message", LVL_MESSAGE},
            {"warning", LVL_WARNING},
            {"error", LVL_ERROR}
    };
}
