#pragma once

#include "tidysq/ProtoSequence.h"

#define IS_PACKED true
#define IS_UNPACKED false
#define IS_CONST true
#define IS_NONCONST false

namespace tidysq {
    namespace internal {
        template<typename INTERNAL, typename PROTO, bool PACKED, bool CONST>
        struct AccessTypeToElementMapper {
            inline static auto map(
                    typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, CONST>::ProtoOrNotSqContentAccessType element_access
                    ) -> typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, CONST>::ProtoOrNotSequenceType {
                return element_access;
            }
        };


        template<typename PROTO, bool CONST>
        struct AccessTypeToElementMapper<RCPP_IT, PROTO, IS_UNPACKED, CONST> {
            inline static auto map(
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_UNPACKED, CONST>::ProtoOrNotSqContentAccessType element_access
            ) -> typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_UNPACKED, CONST>::ProtoOrNotSequenceType {
                return tidysq::ProtoSequence<RCPP_IT, PROTO>(
                        typename util::TypeBinder<RCPP_IT, PROTO>::ProtoSequenceContentStorageType(element_access));
            }
        };

        template<typename PROTO, bool CONST>
        struct AccessTypeToElementMapper<RCPP_IT, PROTO, IS_PACKED, CONST> {
            inline static auto map(
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_PACKED, CONST>::ProtoOrNotSqContentAccessType element_access
            ) -> typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_PACKED, CONST>::ProtoOrNotSequenceType {
                Rcpp::RawVector ret(element_access);
                return tidysq::Sequence<RCPP_IT>(ret, Rcpp::IntegerVector(ret.attr("original_length"))[0]);
            }
        };

        template<typename INTERNAL, typename PROTO, bool PACKED>
        struct AccessTypeAssigner {
            inline static void assign(
                    typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, IS_NONCONST>::ProtoOrNotSqContentAccessType &element_access,
                    typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, IS_NONCONST>::ProtoOrNotSequenceType element
                    ) {
                element_access = element;
            }
        };

        template<typename PROTO>
        struct AccessTypeAssigner<RCPP_IT, PROTO, IS_PACKED> {
            inline static void assign(
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_PACKED, IS_NONCONST>::ProtoOrNotSqContentAccessType &element_access,
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_PACKED, IS_NONCONST>::ProtoOrNotSequenceType element
            ) {
                Rcpp::RawVector content = element.content();
                content.attr("original_length") = element.original_length();
                element_access = content;
            }
        };

        template<typename PROTO>
        struct AccessTypeAssigner<RCPP_IT, PROTO, IS_UNPACKED> {
            inline static void assign(
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_UNPACKED, IS_NONCONST>::ProtoOrNotSqContentAccessType &element_access,
                    typename util::UniversalTypeBinder<RCPP_IT, PROTO, IS_UNPACKED, IS_NONCONST>::ProtoOrNotSequenceType element
            ) {
                element_access = element.content();
            }
        };

        template<typename INTERNAL, typename PROTO, bool PACKED, bool CONST>
        class BasicElementProxy {
            typedef typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, CONST>::ProtoOrNotSqContentAccessType ElementAccessType;
            typedef typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, CONST>::ProtoOrNotSequenceType ElementType;
            typedef typename util::UniversalTypeBinder<INTERNAL, PROTO, PACKED, CONST>::ProtoOrNotSqType ContainerType;

            ElementAccessType contained_element_access_;

            inline explicit BasicElementProxy(ElementAccessType contained_sequence) :
                    contained_element_access_(contained_sequence) {} ;

        public:
            inline operator ElementType() const {
                return AccessTypeToElementMapper<INTERNAL, PROTO, PACKED, CONST>::map(contained_element_access_);
            }

            [[nodiscard]] inline ElementType get() const {
                return this->operator ElementType();
            }

            [[nodiscard]] inline bool operator==(const BasicElementProxy<INTERNAL, PROTO, PACKED, CONST> &other) const {
                return this->get() == other.get();
            }

            [[nodiscard]] inline bool operator==(const BasicElementProxy<INTERNAL, PROTO, PACKED, !CONST> &other) const {
                return this->get() == other.get();
            }

            [[nodiscard]] inline bool operator!=(const BasicElementProxy<INTERNAL, PROTO, PACKED, CONST> &other) const {
                return !operator==(other);
            }

            [[nodiscard]] inline bool operator!=(const BasicElementProxy<INTERNAL, PROTO, PACKED, !CONST> &other) const {
                return !operator==(other);
            }

            template<bool ENABLED = !CONST>
            inline std::enable_if_t<ENABLED, BasicElementProxy&> operator=(const ElementType &element) {
                AccessTypeAssigner<INTERNAL, PROTO, PACKED>::assign(contained_element_access_, element);
                return *this;
            }

            template<bool ENABLED = !CONST>
            inline std::enable_if_t<!ENABLED, BasicElementProxy&> operator=(const ElementType &element) {}

            friend ContainerType;
        };
    }

    template<typename INTERNAL, typename PROTO>
    using ProtoSequenceProxy = internal::BasicElementProxy<INTERNAL, PROTO,  IS_UNPACKED, IS_NONCONST>;

    template <typename INTERNAL>
    using SequenceProxy = internal::BasicElementProxy<INTERNAL, RAWS_PT,  IS_PACKED, IS_NONCONST>;

    template <typename INTERNAL, typename PROTO>
    using ProtoSequenceConstProxy = internal::BasicElementProxy<INTERNAL, PROTO,  IS_UNPACKED, IS_CONST>;

    template <typename INTERNAL>
    using SequenceConstProxy = internal::BasicElementProxy<INTERNAL, RAWS_PT,  IS_PACKED, IS_CONST>;
}

#undef IS_PACKED
#undef IS_UNPACKED
#undef IS_CONST
#undef IS_NONCONST