#pragma once

#include "tidysq/ProtoSequence.h"

namespace tidysq {
    namespace internal {
        template<typename INTERNAL, typename PROTO, bool PACKED, bool CONST>
        class AbstractGeneralSequenceProxy {
        protected:
            typedef typename AccessTypeMapper<INTERNAL, PROTO, PACKED, CONST>::AccessType AccessType;
            typedef typename AccessTypeMapper<INTERNAL, PROTO, PACKED, CONST>::SequenceType SequenceType;
            typedef typename AccessTypeMapper<INTERNAL, PROTO, PACKED, CONST>::SqType SqType;

            AccessType contained_sequence_;

            inline explicit AbstractGeneralSequenceProxy(AccessType contained_sequence) :
                    contained_sequence_(contained_sequence) {} ;

        public:
            inline operator SequenceType() const {
                return contained_sequence_;
            }

            [[nodiscard]] inline SequenceType get() const {
                return (SequenceType) *this;
            }

            [[nodiscard]] inline bool operator==(const AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, CONST> &other) const {
                return this->get() == other.get();
            }

            [[nodiscard]] inline bool operator==(const AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, !CONST> &other) const {
                return this->get() == other.get();
            }

            [[nodiscard]] inline bool operator!=(const AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, CONST> &other) const {
                return !operator==(other);
            }

            [[nodiscard]] inline bool operator!=(const AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, !CONST> &other) const {
                return !operator==(other);
            }

            friend SqType;
        };

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, false, false>::operator ProtoSequence<RCPP_IT, RAWS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, RAWS_PT>(Rcpp::RawVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, false, true>::operator ProtoSequence<RCPP_IT, RAWS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, RAWS_PT>(Rcpp::RawVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, INTS_PT, false, false>::operator ProtoSequence<RCPP_IT, INTS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, INTS_PT>(Rcpp::IntegerVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, INTS_PT, false, true>::operator ProtoSequence<RCPP_IT, INTS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, INTS_PT>(Rcpp::IntegerVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, STRINGS_PT, false, false>::operator ProtoSequence<RCPP_IT, STRINGS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, STRINGS_PT>(Rcpp::StringVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, STRINGS_PT, false, true>::operator ProtoSequence<RCPP_IT, STRINGS_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, STRINGS_PT>(Rcpp::StringVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, STRING_PT, false, false>::operator ProtoSequence<RCPP_IT, STRING_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, STRING_PT>(std::string(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, STRING_PT, false, true>::operator ProtoSequence<RCPP_IT, STRING_PT>() const {
            return tidysq::ProtoSequence<RCPP_IT, STRING_PT>(std::string(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, true, false>::operator Sequence<RCPP_IT>() const {
            Rcpp::RawVector ret(contained_sequence_);
            return tidysq::Sequence<RCPP_IT>(ret, Rcpp::IntegerVector(ret.attr("original_length"))[0]);
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, true, true>::operator Sequence<RCPP_IT>() const {
            Rcpp::RawVector ret(contained_sequence_);
            return tidysq::Sequence<RCPP_IT>(ret, Rcpp::IntegerVector(ret.attr("original_length"))[0]);
        }
    }

    template<typename INTERNAL, typename PROTO, bool PACKED>
    class GeneralSequenceConstProxy : public internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, true> {
        using internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, true>::AbstractGeneralSequenceProxy;

        friend SqType;
    };

    template<typename INTERNAL, typename PROTO, bool PACKED>
    class GeneralSequenceProxy : public internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, false> {
        using internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, false>::AbstractGeneralSequenceProxy;

        friend SqType;

    public:
        inline GeneralSequenceProxy& operator=(const typename internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, false>::SequenceType &general_sequence) {
            internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, false>::contained_sequence_ = general_sequence;
            return *this;
        }
    };

    template<>
    inline GeneralSequenceProxy<RCPP_IT, RAWS_PT, true> &  GeneralSequenceProxy<RCPP_IT, RAWS_PT, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, true, false>::SequenceType &general_sequence) {
        Rcpp::RawVector content = general_sequence.content();
        content.attr("original_length") = general_sequence.original_length();
        contained_sequence_ = content;
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, RAWS_PT, false> &  GeneralSequenceProxy<RCPP_IT, RAWS_PT, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, RAWS_PT, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, INTS_PT, true> &  GeneralSequenceProxy<RCPP_IT, INTS_PT, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, INTS_PT, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, INTS_PT, false> &  GeneralSequenceProxy<RCPP_IT, INTS_PT, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, INTS_PT, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, STRINGS_PT, true> &  GeneralSequenceProxy<RCPP_IT, STRINGS_PT, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, STRINGS_PT, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, STRINGS_PT, false> &  GeneralSequenceProxy<RCPP_IT, STRINGS_PT, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, STRINGS_PT, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, STRING_PT, true> &  GeneralSequenceProxy<RCPP_IT, STRING_PT, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, STRING_PT, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP_IT, STRING_PT, false> &  GeneralSequenceProxy<RCPP_IT, STRING_PT, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP_IT, STRING_PT, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template <typename INTERNAL, typename PROTO>
    using ProtoSequenceProxy = GeneralSequenceProxy<INTERNAL, PROTO, false>;

    template <typename INTERNAL>
    using SequenceProxy = GeneralSequenceProxy<INTERNAL, RAWS_PT, true>;

    template <typename INTERNAL, typename PROTO>
    using ProtoSequenceConstProxy = GeneralSequenceConstProxy<INTERNAL, PROTO, false>;

    template <typename INTERNAL>
    using SequenceConstProxy = GeneralSequenceConstProxy<INTERNAL, RAWS_PT, true>;

}