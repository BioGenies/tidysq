#pragma once

#include "tidysq/types/general.h"

namespace tidysq {
    template<InternalType>
    class Sequence;
}

#include <RcppCommon.h>

namespace Rcpp {
    template<>
    SEXP wrap(const tidysq::Sequence<tidysq::RCPP> &);

    namespace traits {
        template<>
        class Exporter<tidysq::Sequence<tidysq::RCPP>>;
    }
}

#include "tidysq/types/TypeMapper.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/types/SequenceIterator.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sequence {
        typename InternalTypeMapper<INTERNAL>::SequenceContentType content_;
        LenSq original_length_;
    public:
        typedef typename InternalTypeMapper<INTERNAL>::SequenceContentType ContentType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceElementType ElementType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceAccessType AccessType;
        typedef typename InternalTypeMapper<INTERNAL>::SequenceConstAccessType ConstAccessType;

        Sequence(const ContentType &content, const LenSq originalLength) :
                content_(content),
                original_length_(originalLength) {};

        Sequence(const LenSq contentLength, const LenSq originalLength) :
                Sequence(ContentType(contentLength), originalLength) {};

        Sequence() :
                Sequence(0, 0) {};

        Sequence(const Sequence &other) = default;

        Sequence(Sequence &&other) = default;

        Sequence& operator=(const Sequence &other) noexcept = default;

        Sequence& operator=(Sequence &&other) noexcept = default;

        inline AccessType operator[](const LenSq index) {
            return content_[index];
        }

        inline ConstAccessType operator[](const LenSq index) const {
             return content_[index];
        }
        
        SequenceIterator<INTERNAL> begin(const Alphabet& alph) const {
            return SequenceIterator<INTERNAL>(content_, originalLength_, alph);
        }
        
        SequenceIterator<INTERNAL> end(const Alphabet& alph) const {
            return SequenceIterator<INTERNAL>(content_, originalLength_, alph, originalLength_);
        }

        [[nodiscard]] inline LenSq originalLength() const {
            return original_length_;
        }

        [[nodiscard]] inline LenSq length() const {
            return content_.size();
        }

        [[nodiscard]] inline ContentType content() const {
            return content_;
        }

        void trim(const LenSq packed_length, const Alphabet &alphabet) {
            content_.erase(content_.begin() + internal::calculate_packed_internal_length(packed_length, alphabet), content_.end());
            original_length_ = packed_length;
        }

    };
}

namespace Rcpp {
    template<>
    inline SEXP wrap(const tidysq::Sequence<tidysq::RCPP> &obj) {
        Rcpp::RawVector ret = obj.content();
        ret.attr("original_length") = obj.originalLength();
        return ret;
    }

    namespace traits {
        template<>
        class Exporter<tidysq::Sequence<tidysq::RCPP>> {
            typedef typename tidysq::Sequence<tidysq::RCPP> OUT ;
            Rcpp::RawVector vec_;

        public:
            explicit Exporter(SEXP x) :
                    vec_(x) {
                if (TYPEOF(x) != RAWSXP)
                    throw std::invalid_argument("Wrong R type!");
                if (!vec_.hasAttribute("original_length"))
                    throw std::invalid_argument("No 'original_length' attribute!");
            }

            OUT get() {
                OUT ret(OUT(vec_, vec_.attr("original_length")));
                return ret;
            }
        };
    }
}
