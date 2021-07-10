#pragma once

#include <string>
#include <unordered_map>
#include <list>

#include <Rcpp.h>

namespace tidysq {
    typedef R_xlen_t                LenSq;
    typedef unsigned long int       u_LenSq;          //max vector size for comparisons
    typedef unsigned char           ElementPacked;
    typedef unsigned char           ElementRaws;
    typedef unsigned short int      ElementInts;
    typedef std::string             ElementStrings;
    typedef char                    ElementStringSimple;
    typedef std::string             ElementStringMultichar;
    typedef unsigned short int      AlphSize;
    typedef unsigned short int      LetterValue;
    typedef std::string             Letter;
    typedef char                    SimpleLetter;

    template<typename INTERNAL, typename PROTO>
    class ProtoSequence;

    template<typename INTERNAL, typename PROTO>
    class ProtoSq;

    template<typename INTERNAL>
    class Sequence;

    template<typename INTERNAL>
    class Sq;

    struct InternalType {
        // structure that extends InternalType should contain following typedefs:
        // SequenceElementType          - type of single packed element, intentionally it should be ElementPacked;
        // SequenceContentStorageType   - type of sequence content storage type;
        // SequenceContentAccessType           - decltype of SequenceContentStorageType::operator[]
        // SequenceContentConstAccessType      - decltype of SequenceContentStorageType::operator[] const
        // SqContentStorageType         - type of sq content storage type
    };

    struct STD_IT : public InternalType {
        typedef Sequence<STD_IT>                    SequenceType;
        typedef ElementPacked                       SequenceElementType;
        typedef std::vector<ElementPacked>          SequenceContentStorageType;
        typedef ElementPacked &                     SequenceContentAccessType;
        typedef const ElementPacked &               SequenceContentConstAccessType;

        typedef Sq<STD_IT>                          SqType;
        typedef std::vector<Sequence<STD_IT>>       SqContentStorageType;
        typedef Sequence<STD_IT> &                  SqContentAccessType;
        typedef const Sequence<STD_IT> &            SqContentConstAccessType;
    };

    struct RCPP_IT : public InternalType {
        typedef Sequence<RCPP_IT>                   SequenceType;
        typedef ElementPacked                       SequenceElementType;
        typedef Rcpp::RawVector                     SequenceContentStorageType;
        typedef Rcpp::RawVector::Proxy              SequenceContentAccessType;
        typedef Rcpp::RawVector::const_Proxy        SequenceContentConstAccessType;

        typedef Sq<RCPP_IT>                         SqType;
        typedef Rcpp::List                          SqContentStorageType;
        typedef Rcpp::List::Proxy                   SqContentAccessType;
        typedef Rcpp::List::const_Proxy             SqContentConstAccessType;
    };

    struct ProtoType {};

    struct RAWS_PT : public ProtoType {
        typedef ElementRaws                                 ProtoSequenceElementType;
    };

    struct INTS_PT : public ProtoType {
        typedef ElementInts                                 ProtoSequenceElementType;
    };

    struct STRINGS_PT : public ProtoType {
        typedef ElementStrings                              ProtoSequenceElementType;
    };

    struct STRING_PT : public ProtoType {
        typedef ElementStringSimple                         ProtoSequenceElementType;
    };

#define COMMON_STRING_PROTO_TYPE_TYPEDEFS_ \
    typedef std::string                                                         ProtoSequenceContentStorageType; \
    typedef ElementStringSimple &                                               ProtoSequenceContentAccessType; \
    typedef const ElementStringSimple &                                         ProtoSequenceContentConstAccessType;

#define COMMON_RCPP_INTERNAL_TYPE_TYPEDEFS_ \
    typedef Rcpp::List                                                          ProtoSqListConstructorType; \
    typedef Rcpp::List                                                          ProtoSqContentStorageType; \
    typedef Rcpp::List::Proxy                                                   ProtoSqContentAccessType; \
    typedef Rcpp::List::const_Proxy                                             ProtoSqContentConstAccessType;

#define COMMON_RCPP_INTERNAL_TYPE_PROTO_SEQUENCE_TYPEDEFS_(VECTOR_TYPE) \
    typedef VECTOR_TYPE                                                         ProtoSequenceContentStorageType; \
    typedef VECTOR_TYPE::Proxy                                                  ProtoSequenceContentAccessType; \
    typedef VECTOR_TYPE::const_Proxy                                            ProtoSequenceContentConstAccessType;

    namespace util {
        template<typename INTERNAL, typename PROTO>
        struct TypeBinder;

        template<typename PROTO>
        struct TypeBinder<STD_IT, PROTO> {
            typedef ProtoSequence<STD_IT, PROTO>                                        ProtoSequenceType;
            typedef std::vector<typename PROTO::ProtoSequenceElementType>               ProtoSequenceContentStorageType;
            typedef typename PROTO::ProtoSequenceElementType &                          ProtoSequenceContentAccessType;
            typedef const typename PROTO::ProtoSequenceElementType &                    ProtoSequenceContentConstAccessType;

            typedef ProtoSq<STD_IT, PROTO>                                              ProtoSqType;
            typedef std::vector<std::vector<typename PROTO::ProtoSequenceElementType>>  ProtoSqListConstructorType;
            typedef std::vector<ProtoSequence<STD_IT, PROTO>>                           ProtoSqContentStorageType;
            typedef ProtoSequence<STD_IT, PROTO> &                                      ProtoSqContentAccessType;
            typedef const ProtoSequence<STD_IT, PROTO> &                                ProtoSqContentConstAccessType;
        };

        template<>
        struct TypeBinder<STD_IT, STRING_PT> {
            typedef ProtoSequence<STD_IT, STRING_PT>                                    ProtoSequenceType;
            COMMON_STRING_PROTO_TYPE_TYPEDEFS_

            typedef ProtoSq<STD_IT, STRING_PT>                                          ProtoSqType;
            typedef std::vector<std::string>                                            ProtoSqListConstructorType;
            typedef std::vector<ProtoSequence<STD_IT, STRING_PT>>                       ProtoSqContentStorageType;
            typedef ProtoSequence<STD_IT, STRING_PT> &                                  ProtoSqContentAccessType;
            typedef const ProtoSequence<STD_IT, STRING_PT> &                            ProtoSqContentConstAccessType;
        };

        template<>
        struct TypeBinder<RCPP_IT, RAWS_PT> {
            typedef ProtoSequence<RCPP_IT, RAWS_PT>                                     ProtoSequenceType;
            COMMON_RCPP_INTERNAL_TYPE_PROTO_SEQUENCE_TYPEDEFS_(Rcpp::RawVector)

            typedef ProtoSq<RCPP_IT, RAWS_PT>                                           ProtoSqType;
            COMMON_RCPP_INTERNAL_TYPE_TYPEDEFS_
        };

        template<>
        struct TypeBinder<RCPP_IT, INTS_PT> {
            typedef ProtoSequence<RCPP_IT, INTS_PT>                                     ProtoSequenceType;
            COMMON_RCPP_INTERNAL_TYPE_PROTO_SEQUENCE_TYPEDEFS_(Rcpp::IntegerVector)

            typedef ProtoSq<RCPP_IT, INTS_PT>                                           ProtoSqType;
            COMMON_RCPP_INTERNAL_TYPE_TYPEDEFS_
        };

        template<>
        struct TypeBinder<RCPP_IT, STRINGS_PT> {
            typedef ProtoSequence<RCPP_IT, STRINGS_PT>                                  ProtoSequenceType;
            COMMON_RCPP_INTERNAL_TYPE_PROTO_SEQUENCE_TYPEDEFS_(Rcpp::StringVector)

            typedef ProtoSq<RCPP_IT, STRINGS_PT>                                        ProtoSqType;
            COMMON_RCPP_INTERNAL_TYPE_TYPEDEFS_
        };

        template<>
        struct TypeBinder<RCPP_IT, STRING_PT> {
            typedef ProtoSequence<RCPP_IT, STRING_PT>                                   ProtoSequenceType;
            COMMON_STRING_PROTO_TYPE_TYPEDEFS_

            typedef ProtoSq<RCPP_IT, STRING_PT>                                         ProtoSqType;
            typedef Rcpp::StringVector                                                  ProtoSqListConstructorType;
            typedef Rcpp::StringVector                                                  ProtoSqContentStorageType;
            typedef Rcpp::StringVector::Proxy                                           ProtoSqContentAccessType;
            typedef Rcpp::StringVector::const_Proxy                                     ProtoSqContentConstAccessType;
        };

        template<typename INTERNAL, typename PROTO, bool PACKED, bool CONST>
        struct UniversalTypeBinder {
            typedef std::conditional_t<PACKED,
                    typename INTERNAL::SequenceType,
                    typename TypeBinder<INTERNAL, PROTO>::ProtoSequenceType>                        ProtoOrNotSequenceType;
            typedef std::conditional_t<PACKED,
                    typename INTERNAL::SequenceContentStorageType,
                    typename TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentStorageType>          ProtoOrNotSequenceContentStorageType;
            typedef std::conditional_t<PACKED,
                    std::conditional_t<CONST,
                            typename INTERNAL::SequenceContentConstAccessType,
                            typename INTERNAL::SequenceContentAccessType>,
                    std::conditional_t<CONST,
                            typename TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentConstAccessType,
                            typename TypeBinder<INTERNAL, PROTO>::ProtoSequenceContentAccessType>>  ProtoOrNotSequenceContentAccessType;

            typedef std::conditional_t<PACKED,
                    typename INTERNAL::SqType,
                    typename TypeBinder<INTERNAL, PROTO>::ProtoSqType>                              ProtoOrNotSqType;
            typedef std::conditional_t<PACKED,
                    typename INTERNAL::SqContentStorageType,
                    typename TypeBinder<INTERNAL, PROTO>::ProtoSqContentStorageType>                ProtoOrNotSqContentStorageType;
            typedef std::conditional_t<PACKED,
                    std::conditional_t<CONST,
                            typename INTERNAL::SqContentConstAccessType,
                            typename INTERNAL::SqContentAccessType>,
                    std::conditional_t<CONST,
                            typename TypeBinder<INTERNAL, PROTO>::ProtoSqContentConstAccessType,
                            typename TypeBinder<INTERNAL, PROTO>::ProtoSqContentAccessType>>        ProtoOrNotSqContentAccessType;
        };
    }

#undef COMMON_STRING_PROTO_TYPE_TYPEDEFS_
#undef COMMON_RCPP_INTERNAL_TYPE_TYPEDEFS_
#undef COMMON_RCPP_INTERNAL_TYPE_PROTO_SEQUENCE_TYPEDEFS_

    enum SqType {
        AMI_EXT,
        AMI_BSC,
        DNA_EXT,
        DNA_BSC,
        RNA_EXT,
        RNA_BSC,
        UNT,
        ATP,
        ENC
    };

    namespace internal {
        typedef std::unordered_map<LetterValue, const LetterValue>      ComplementTable;
        typedef std::unordered_map<LetterValue,
                const std::unordered_map<LetterValue,
                        const std::unordered_map<LetterValue,
                                const LetterValue>>>                    CodonTable;
        typedef std::unordered_map<Letter, std::list<Letter>>           AmbiguousDict;
    }
}
