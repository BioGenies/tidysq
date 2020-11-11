#pragma once

#include "tidysq/ProtoSequence.h"

namespace tidysq {
    namespace internal {
        template<InternalType INTERNAL, ProtoType PROTO, bool PACKED, bool CONST>
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
        inline AbstractGeneralSequenceProxy<RCPP, RAWS, false, false>::operator ProtoSequence<RCPP, RAWS>() const {
            return tidysq::ProtoSequence<RCPP, RAWS>(Rcpp::RawVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, RAWS, false, true>::operator ProtoSequence<RCPP, RAWS>() const {
            return tidysq::ProtoSequence<RCPP, RAWS>(Rcpp::RawVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, INTS, false, false>::operator ProtoSequence<RCPP, INTS>() const {
            return tidysq::ProtoSequence<RCPP, INTS>(Rcpp::IntegerVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, INTS, false, true>::operator ProtoSequence<RCPP, INTS>() const {
            return tidysq::ProtoSequence<RCPP, INTS>(Rcpp::IntegerVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, STRINGS, false, false>::operator ProtoSequence<RCPP, STRINGS>() const {
            return tidysq::ProtoSequence<RCPP, STRINGS>(Rcpp::StringVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, STRINGS, false, true>::operator ProtoSequence<RCPP, STRINGS>() const {
            return tidysq::ProtoSequence<RCPP, STRINGS>(Rcpp::StringVector(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, STRING, false, false>::operator ProtoSequence<RCPP, STRING>() const {
            return tidysq::ProtoSequence<RCPP, STRING>(std::string(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, STRING, false, true>::operator ProtoSequence<RCPP, STRING>() const {
            return tidysq::ProtoSequence<RCPP, STRING>(std::string(contained_sequence_));
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, RAWS, true, false>::operator Sequence<RCPP>() const {
            Rcpp::RawVector ret(contained_sequence_);
            return tidysq::Sequence<RCPP>(ret, Rcpp::IntegerVector(ret.attr("original_length"))[0]);
        }

        template<>
        inline AbstractGeneralSequenceProxy<RCPP, RAWS, true, true>::operator Sequence<RCPP>() const {
            Rcpp::RawVector ret(contained_sequence_);
            return tidysq::Sequence<RCPP>(ret, Rcpp::IntegerVector(ret.attr("original_length"))[0]);
        }
    }

    template<InternalType INTERNAL, ProtoType PROTO, bool PACKED>
    class GeneralSequenceConstProxy : public internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, true> {
        using internal::AbstractGeneralSequenceProxy<INTERNAL, PROTO, PACKED, true>::AbstractGeneralSequenceProxy;

        friend SqType;
    };

    template<InternalType INTERNAL, ProtoType PROTO, bool PACKED>
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
    inline GeneralSequenceProxy<RCPP, RAWS, true> &  GeneralSequenceProxy<RCPP, RAWS, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, RAWS, true, false>::SequenceType &general_sequence) {
        Rcpp::RawVector content = general_sequence.content();
        content.attr("original_length") = general_sequence.originalLength();
        contained_sequence_ = content;
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, RAWS, false> &  GeneralSequenceProxy<RCPP, RAWS, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, RAWS, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, INTS, true> &  GeneralSequenceProxy<RCPP, INTS, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, INTS, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, INTS, false> &  GeneralSequenceProxy<RCPP, INTS, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, INTS, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, STRINGS, true> &  GeneralSequenceProxy<RCPP, STRINGS, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, STRINGS, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, STRINGS, false> &  GeneralSequenceProxy<RCPP, STRINGS, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, STRINGS, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, STRING, true> &  GeneralSequenceProxy<RCPP, STRING, true>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, STRING, true, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template<>
    inline GeneralSequenceProxy<RCPP, STRING, false> &  GeneralSequenceProxy<RCPP, STRING, false>::operator=(
            const typename internal::AbstractGeneralSequenceProxy<RCPP, STRING, false, false>::SequenceType &proto) {
        contained_sequence_ = proto.content();
        return *this;
    }

    template <InternalType INTERNAL, ProtoType PROTO>
    using ProtoSequenceProxy = GeneralSequenceProxy<INTERNAL, PROTO, false>;

    template <InternalType INTERNAL>
    using SequenceProxy = GeneralSequenceProxy<INTERNAL, RAWS, true>;

    template <InternalType INTERNAL, ProtoType PROTO>
    using ProtoSequenceConstProxy = GeneralSequenceConstProxy<INTERNAL, PROTO, false>;

    template <InternalType INTERNAL>
    using SequenceConstProxy = GeneralSequenceConstProxy<INTERNAL, RAWS, true>;

}