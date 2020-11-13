#pragma once

#include <vector>
#include <Rcpp.h>

#include "tidysq/tidysq-typedefs.h"

namespace tidysq {
    template<typename INTERNAL>
    class Sequence;

    template<typename INTERNAL, typename PROTO>
    class ProtoSequence;

    template<typename INTERNAL>
    class Sq;

    template<typename INTERNAL, typename PROTO>
    class ProtoSq;

    template<typename INTERNAL>
    struct InternalTypeMapper;

    template<>
    struct InternalTypeMapper<STD_IT> {
        typedef ElementPacked                                                      SequenceElementType;
        typedef std::vector<ElementPacked>                                         SequenceContentType;
        typedef         SequenceElementType &                                   SequenceAccessType;
        typedef const   SequenceElementType &                                   SequenceConstAccessType;

        typedef std::vector<Sequence<STD_IT>>                                      SqContentType;
        typedef Sequence<STD_IT>                                                   SqElementType;
    };

    template<>
    struct InternalTypeMapper<RCPP_IT> {
        typedef ElementPacked                                                      SequenceElementType;
        typedef Rcpp::RawVector                                                 SequenceContentType;
        typedef Rcpp::RawVector::Proxy                                          SequenceAccessType;
        typedef Rcpp::RawVector::const_Proxy                                    SequenceConstAccessType;

        typedef Rcpp::List                                                      SqContentType;
        typedef Sequence<RCPP_IT>                                                  SqElementType;
    };

    template<typename PROTO>
    struct ProtoTypeMapper;

    template<>
    struct ProtoTypeMapper<RAWS_PT> {
        typedef ElementRaws                                                     ProtoSequenceElementType;
        typedef std::vector<ElementRaws>                                        ProtoSequenceStdContentType;
        typedef Rcpp::RawVector                                                 ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<INTS_PT> {
        typedef ElementInts                                                     ProtoSequenceElementType;
        typedef std::vector<ElementInts>                                        ProtoSequenceStdContentType;
        typedef Rcpp::IntegerVector                                             ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRINGS_PT> {
        typedef ElementStrings                                                  ProtoSequenceElementType;
        typedef std::vector<ElementStrings>                                     ProtoSequenceStdContentType;
        typedef Rcpp::StringVector                                              ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRING_PT> {
        typedef ElementStringSimple                                             ProtoSequenceElementType;
        typedef std::string                                                     ProtoSequenceStdContentType;
        typedef std::string                                                     ProtoSequenceRcppContentType;
    };

    template<typename INTERNAL, typename PROTO>
    struct TypeMapper;
    
    template<typename PROTO>
    struct TypeMapper<STD_IT, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceStdContentType    ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef ProtoSequence<STD_IT, PROTO>                                       ProtoSqElementType;
        typedef std::vector<ProtoSqElementType>                                 ProtoSqContentType;

    };

    template<typename PROTO>
    struct TypeMapper<RCPP_IT, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceRcppContentType   ProtoSequenceContentType;
        typedef typename ProtoSequenceContentType::Proxy                        ProtoSequenceAccessType;
        typedef typename ProtoSequenceContentType::const_Proxy                  ProtoSequenceConstAccessType;

        typedef Rcpp::List                                                      ProtoSqContentType;
        typedef ProtoSequence<RCPP_IT, PROTO>                                      ProtoSqElementType;
    };

    template<>
    struct TypeMapper<RCPP_IT, STRING_PT> {
        typedef typename ProtoTypeMapper<STRING_PT>::ProtoSequenceElementType      ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<STRING_PT>::ProtoSequenceRcppContentType  ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef Rcpp::StringVector                                              ProtoSqContentType;
        typedef ProtoSequence<RCPP_IT, STRING_PT>                                     ProtoSqElementType;
    };


    template<typename INTERNAL, typename PROTO, bool PACKED, bool CONST>
    struct AccessTypeMapper;

    template<typename PROTO>
    struct AccessTypeMapper<STD_IT, PROTO, false, true> {
        typedef const typename TypeMapper<STD_IT, PROTO>::ProtoSqElementType &   AccessType;
        typedef ProtoSequence<STD_IT, PROTO>                                           SequenceType;
        typedef ProtoSq<STD_IT, PROTO>                                                 SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<STD_IT, PROTO, false, false> {
        typedef       typename TypeMapper<STD_IT, PROTO>::ProtoSqElementType &   AccessType;
        typedef ProtoSequence<STD_IT, PROTO>                                           SequenceType;
        typedef ProtoSq<STD_IT, PROTO>                                                 SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<STD_IT, PROTO, true, true> {
        typedef const typename InternalTypeMapper<STD_IT>::SqElementType &       AccessType;
        typedef Sequence<STD_IT>                                                       SequenceType;
        typedef Sq<STD_IT>                                                              SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<STD_IT, PROTO, true, false> {
        typedef       typename InternalTypeMapper<STD_IT>::SqElementType &       AccessType;
        typedef Sequence<STD_IT>                                                       SequenceType;
        typedef Sq<STD_IT>                                                              SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<RCPP_IT, PROTO, false, true> {
        typedef Rcpp::List::const_Proxy                                                   AccessType;
        typedef ProtoSequence<RCPP_IT, PROTO>                                                SequenceType;
        typedef ProtoSq<RCPP_IT, PROTO>                                                     SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<RCPP_IT, PROTO, false, false> {
        typedef Rcpp::List::Proxy         AccessType;
        typedef ProtoSequence<RCPP_IT, PROTO>                                                SequenceType;
        typedef ProtoSq<RCPP_IT, PROTO>                                                     SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<RCPP_IT, PROTO, true, true> {
        typedef Rcpp::List::const_Proxy        AccessType;
        typedef Sequence<RCPP_IT>                                                       SequenceType;
        typedef Sq<RCPP_IT>                                                              SqType;
    };

    template<typename PROTO>
    struct AccessTypeMapper<RCPP_IT, PROTO, true, false> {
        typedef Rcpp::List::Proxy        AccessType;
        typedef Sequence<RCPP_IT>                                                       SequenceType;
        typedef Sq<RCPP_IT>                                                              SqType;
    };

//    template<>
//    struct AccessTypeMapper<RCPP_IT, STRING_PT, true, true> {
//        typedef const typename InternalTypeMapper<RCPP_IT>::SequenceContentType &        AccessType;
//        typedef Sequence<RCPP_IT>                                                       SequenceType;
//    };
//
//    template<>
//    struct AccessTypeMapper<RCPP_IT, STRING_PT, true, false> {
//        typedef      typename InternalTypeMapper<RCPP_IT>::SequenceContentType &        AccessType;
//        typedef Sequence<RCPP_IT>                                                       SequenceType;
//    };

    template<>
    struct AccessTypeMapper<RCPP_IT, STRING_PT, false, true> {
        typedef typename TypeMapper<RCPP_IT, STRING_PT>::ProtoSqContentType::const_Proxy  AccessType;
        typedef ProtoSequence<RCPP_IT, STRING_PT>                                         SequenceType;
        typedef ProtoSq<RCPP_IT, STRING_PT>                                                SqType;
    };

    template<>
    struct AccessTypeMapper<RCPP_IT, STRING_PT, false, false> {
        typedef typename TypeMapper<RCPP_IT, STRING_PT>::ProtoSqContentType::Proxy     AccessType;
        typedef ProtoSequence<RCPP_IT, STRING_PT>                                         SequenceType;
        typedef ProtoSq<RCPP_IT, STRING_PT>                                                SqType;
    };
}

