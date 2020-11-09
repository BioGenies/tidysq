#pragma once

namespace tidysq::internal {
    template<typename TYPE, bool CONST>
    struct ConstWrapper;

    template<typename TYPE>
    struct ConstWrapper<TYPE, false> {
        typedef TYPE   Type;
        typedef TYPE & Reference;
    };

    template<typename TYPE>
    struct ConstWrapper<TYPE, true> {
        typedef const TYPE   Type;
        typedef const TYPE & Reference;
    };
}
