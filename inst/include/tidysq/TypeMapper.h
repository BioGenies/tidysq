#pragma once

#include <vector>
#include <Rcpp.h>

#include "tidysq/tidysq-typedefs.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sequence;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSequence;

    template<InternalType INTERNAL>
    class Sq;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSq;

    template<InternalType INTERNAL>
    struct InternalTypeMapper;

    template<>
    struct InternalTypeMapper<STD> {
        typedef ElementPacked                                                      SequenceElementType;
        typedef std::vector<ElementPacked>                                         SequenceContentType;
        typedef         SequenceElementType &                                   SequenceAccessType;
        typedef const   SequenceElementType &                                   SequenceConstAccessType;

        typedef std::vector<Sequence<STD>>                                      SqContentType;
        typedef Sequence<STD>                                                   SqElementType;
    };

    template<>
    struct InternalTypeMapper<RCPP> {
        typedef ElementPacked                                                      SequenceElementType;
        typedef Rcpp::RawVector                                                 SequenceContentType;
        typedef Rcpp::RawVector::Proxy                                          SequenceAccessType;
        typedef Rcpp::RawVector::const_Proxy                                    SequenceConstAccessType;

        typedef Rcpp::List                                                      SqContentType;
        typedef Sequence<RCPP>                                                  SqElementType;
    };

    template<ProtoType PROTO>
    struct ProtoTypeMapper;

    template<>
    struct ProtoTypeMapper<RAWS> {
        typedef ElementRaws                                                     ProtoSequenceElementType;
        typedef std::vector<ElementRaws>                                        ProtoSequenceStdContentType;
        typedef Rcpp::RawVector                                                 ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<INTS> {
        typedef ElementInts                                                     ProtoSequenceElementType;
        typedef std::vector<ElementInts>                                        ProtoSequenceStdContentType;
        typedef Rcpp::IntegerVector                                             ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRINGS> {
        typedef ElementStrings                                                  ProtoSequenceElementType;
        typedef std::vector<ElementStrings>                                     ProtoSequenceStdContentType;
        typedef Rcpp::StringVector                                              ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRING> {
        typedef ElementStringSimple                                             ProtoSequenceElementType;
        typedef std::string                                                     ProtoSequenceStdContentType;
        typedef std::string                                                     ProtoSequenceRcppContentType;
    };

    template<InternalType INTERNAL, ProtoType PROTO>
    struct TypeMapper;
    
    template<ProtoType PROTO>
    struct TypeMapper<STD, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceStdContentType    ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef ProtoSequence<STD, PROTO>                                       ProtoSqElementType;
        typedef std::vector<ProtoSqElementType>                                 ProtoSqContentType;

    };

    template<ProtoType PROTO>
    struct TypeMapper<RCPP, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceRcppContentType   ProtoSequenceContentType;
        typedef typename ProtoSequenceContentType::Proxy                        ProtoSequenceAccessType;
        typedef typename ProtoSequenceContentType::const_Proxy                  ProtoSequenceConstAccessType;

        typedef Rcpp::List                                                      ProtoSqContentType;
        typedef ProtoSequence<RCPP, PROTO>                                      ProtoSqElementType;
    };

    template<>
    struct TypeMapper<RCPP, STRING> {
        typedef typename ProtoTypeMapper<STRING>::ProtoSequenceElementType      ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<STRING>::ProtoSequenceRcppContentType  ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef Rcpp::StringVector                                              ProtoSqContentType;
        typedef ProtoSequence<RCPP, STRING>                                     ProtoSqElementType;
    };


    template<InternalType INTERNAL, ProtoType PROTO, bool PACKED, bool CONST>
    struct AccessTypeMapper;

    template<ProtoType PROTO>
    struct AccessTypeMapper<STD, PROTO, false, true> {
        typedef const typename TypeMapper<STD, PROTO>::ProtoSqElementType &   AccessType;
        typedef ProtoSequence<STD, PROTO>                                           SequenceType;
        typedef ProtoSq<STD, PROTO>                                                 SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<STD, PROTO, false, false> {
        typedef       typename TypeMapper<STD, PROTO>::ProtoSqElementType &   AccessType;
        typedef ProtoSequence<STD, PROTO>                                           SequenceType;
        typedef ProtoSq<STD, PROTO>                                                 SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<STD, PROTO, true, true> {
        typedef const typename InternalTypeMapper<STD>::SqElementType &       AccessType;
        typedef Sequence<STD>                                                       SequenceType;
        typedef Sq<STD>                                                              SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<STD, PROTO, true, false> {
        typedef       typename InternalTypeMapper<STD>::SqElementType &       AccessType;
        typedef Sequence<STD>                                                       SequenceType;
        typedef Sq<STD>                                                              SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<RCPP, PROTO, false, true> {
        typedef Rcpp::List::const_Proxy                                                   AccessType;
        typedef ProtoSequence<RCPP, PROTO>                                                SequenceType;
        typedef ProtoSq<RCPP, PROTO>                                                     SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<RCPP, PROTO, false, false> {
        typedef Rcpp::List::Proxy         AccessType;
        typedef ProtoSequence<RCPP, PROTO>                                                SequenceType;
        typedef ProtoSq<RCPP, PROTO>                                                     SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<RCPP, PROTO, true, true> {
        typedef Rcpp::List::const_Proxy        AccessType;
        typedef Sequence<RCPP>                                                       SequenceType;
        typedef Sq<RCPP>                                                              SqType;
    };

    template<ProtoType PROTO>
    struct AccessTypeMapper<RCPP, PROTO, true, false> {
        typedef Rcpp::List::Proxy        AccessType;
        typedef Sequence<RCPP>                                                       SequenceType;
        typedef Sq<RCPP>                                                              SqType;
    };

//    template<>
//    struct AccessTypeMapper<RCPP, STRING, true, true> {
//        typedef const typename InternalTypeMapper<RCPP>::SequenceContentType &        AccessType;
//        typedef Sequence<RCPP>                                                       SequenceType;
//    };
//
//    template<>
//    struct AccessTypeMapper<RCPP, STRING, true, false> {
//        typedef      typename InternalTypeMapper<RCPP>::SequenceContentType &        AccessType;
//        typedef Sequence<RCPP>                                                       SequenceType;
//    };

    template<>
    struct AccessTypeMapper<RCPP, STRING, false, true> {
        typedef typename TypeMapper<RCPP, STRING>::ProtoSqContentType::const_Proxy  AccessType;
        typedef ProtoSequence<RCPP, STRING>                                         SequenceType;
        typedef ProtoSq<RCPP, STRING>                                                SqType;
    };

    template<>
    struct AccessTypeMapper<RCPP, STRING, false, false> {
        typedef typename TypeMapper<RCPP, STRING>::ProtoSqContentType::Proxy     AccessType;
        typedef ProtoSequence<RCPP, STRING>                                         SequenceType;
        typedef ProtoSq<RCPP, STRING>                                                SqType;
    };
}

