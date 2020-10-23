#pragma once

#include "tidysq/types/general.h"

namespace tidysq {
    template<InternalType, ProtoType>
    class ProtoSequence;
}

#include <RcppCommon.h>

namespace Rcpp {
    template<>
    SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::RAWS> &);

    template<>
    SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::INTS> &);

    template<>
    SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::STRINGS> &);

    template<>
    SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::STRING> &);
}

#include <utility>

#include "tidysq/types/Alphabet.h"
#include "tidysq/types/TypeMapper.h"
#include "tidysq/ops/internal/util.h"
#include "tidysq/types/ProtoSequenceInputInterpreter.h"

namespace tidysq {
    template<InternalType INTERNAL, ProtoType PROTO, bool SIMPLE>
    class ProtoSequenceInputInterpreter;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSequence {
        typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType content_;
    public:
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceContentType ContentType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceElementType ElementType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceAccessType AccessType;
        typedef typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceConstAccessType ConstAccessType;

        explicit ProtoSequence(const ContentType &content) :
                content_(content) {};

        explicit ProtoSequence(const LenSq length) :
                ProtoSequence(ContentType(length)) {};

        ProtoSequence() :
                ProtoSequence(0) {};

        ProtoSequence(const std::initializer_list<ElementType> &list) :
                content_(ContentType(list)) {};

        ProtoSequence(const ProtoSequence &other) = default;

        ProtoSequence(ProtoSequence &&other)  noexcept = default;

        ProtoSequence& operator=(const ProtoSequence &other) noexcept = default;

        ProtoSequence& operator=(ProtoSequence &&other) noexcept = default;

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        };

        inline ConstAccessType operator[](const LenSq index) const {
            return content_[index];
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline const ContentType &content() const {
            return content_;
        }

        inline bool operator==(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return content_ == other.content_;
        }

        inline bool operator!=(const ProtoSequence<INTERNAL, PROTO> &other) const {
            return !operator==(other);
        }

        inline ProtoSequence<INTERNAL, PROTO>& operator+=(std::string &&letter) {
            throw std::exception();
        }

        template<bool SIMPLE>
        inline ProtoSequenceInputInterpreter<INTERNAL, PROTO, SIMPLE> content_interpreter(const Alphabet &alphabet) const {
            return ProtoSequenceInputInterpreter<INTERNAL, PROTO, SIMPLE>(content_.cbegin(), content_.cend(), alphabet);
        }
    };

    template<>
    inline ProtoSequence<STD, STRING>::ProtoSequence(const LenSq length) :
            ProtoSequence(ContentType(length, ' ')) {};

    template<>
    inline ProtoSequence<RCPP, STRING>::ProtoSequence(const LenSq length) :
            ProtoSequence(ContentType(length, ' ')) {};


    template<>
    inline ProtoSequence<RCPP, STRING>& ProtoSequence<RCPP, STRING>::operator+=(std::string &&letter) {
        content_ += letter;
        return *this;
    }

    template<>
    inline ProtoSequence<STD, STRING>& ProtoSequence<STD, STRING>::operator+=(std::string &&letter) {
        std::string x;
        content_ += letter;
        return *this;
    }
}

namespace Rcpp {
    template<>
    inline SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::RAWS> &obj) {
        return obj.content();
    }

    template<>
    inline SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::INTS> &obj) {
        return obj.content();
    }

    template<>
    inline SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::STRINGS> &obj) {
        return obj.content();
    }

    template<>
    inline SEXP wrap(const tidysq::ProtoSequence<tidysq::RCPP, tidysq::STRING> &obj) {
        return Rcpp::StringVector(obj.content());
    }
}
