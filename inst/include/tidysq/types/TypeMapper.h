#ifndef TIDYSQ_TYPEMAPPER_H
#define TIDYSQ_TYPEMAPPER_H

#include <vector>
#include <Rcpp.h>

#include "general.h"

namespace tidysq {
    template<InternalType INTERNAL>
    class Sequence;

    template<InternalType INTERNAL, ProtoType PROTO>
    class ProtoSequence;

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
        typedef         Sequence<STD> &                                         SqAccessType;
        typedef const   Sequence<STD> &                                         SqConstAccessType;
    };

    template<>
    struct InternalTypeMapper<RCPP> {
        typedef ElementPacked                                                      SequenceElementType;
        typedef Rcpp::RawVector                                                 SequenceContentType;
        typedef Rcpp::RawVector::Proxy                                          SequenceAccessType;
        typedef Rcpp::RawVector::const_Proxy                                    SequenceConstAccessType;

        typedef Rcpp::List                                                      SqContentType;
        typedef Sequence<RCPP>                                                  SqElementType;
        typedef         Rcpp::List::Proxy                                       SqAccessType;
        typedef const   Rcpp::List::const_Proxy                                 SqConstAccessType;
    };

    template<ProtoType PROTO>
    struct ProtoTypeMapper;

    template<>
    struct ProtoTypeMapper<RAWS> {
        typedef ElementRaws                                                        ProtoSequenceElementType;
        typedef Rcpp::RawVector                                                 ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<INTS> {
        typedef ElementInts                                                        ProtoSequenceElementType;
        typedef Rcpp::IntegerVector                                             ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRINGS> {
        typedef ElementStrings                                                     ProtoSequenceElementType;
        typedef Rcpp::StringVector                                              ProtoSequenceRcppContentType;
    };

    template<>
    struct ProtoTypeMapper<STRING> {
        typedef ElementStringSimple                                                      ProtoSequenceElementType;
        typedef std::string                                                     ProtoSequenceRcppContentType;
    };

    template<InternalType INTERNAL, ProtoType PROTO>
    struct TypeMapper;
    
    template<ProtoType PROTO>
    struct TypeMapper<STD, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef std::vector<ProtoSequenceElementType>                           ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef ProtoSequence<STD, RAWS>                                        ProtoSqElementType;
        typedef std::vector<ProtoSqElementType>                                 ProtoSqContentType;
        typedef         ProtoSqElementType &                                    ProtoSqAccessType;
        typedef const   ProtoSqElementType &                                    ProtoSqConstAccessType;

    };

    template<ProtoType PROTO>
    struct TypeMapper<RCPP, PROTO> {
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceElementType       ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<PROTO>::ProtoSequenceRcppContentType   ProtoSequenceContentType;
        typedef typename ProtoSequenceContentType::Proxy                        ProtoSequenceAccessType;
        typedef typename ProtoSequenceContentType::const_Proxy                  ProtoSequenceConstAccessType;

        typedef Rcpp::List                                                      ProtoSqContentType;
        typedef ProtoSequence<RCPP, PROTO>                                      ProtoSqElementType;
        typedef         ProtoSqContentType::Proxy                               ProtoSqAccessType;
        typedef const   ProtoSqContentType::const_Proxy                         ProtoSqConstAccessType;
    };

    template<>
    struct TypeMapper<RCPP, STRING> {
        typedef typename ProtoTypeMapper<STRING>::ProtoSequenceElementType      ProtoSequenceElementType;
        typedef typename ProtoTypeMapper<STRING>::ProtoSequenceRcppContentType  ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &                              ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &                              ProtoSequenceConstAccessType;

        typedef Rcpp::StringVector                                              ProtoSqContentType;
        typedef ProtoSequence<RCPP, STRING>                                     ProtoSqElementType;
        typedef         ProtoSqContentType::Proxy                               ProtoSqAccessType;
        typedef const   ProtoSqContentType::const_Proxy                         ProtoSqConstAccessType;
    };
}

#endif //TIDYSQ_TYPEMAPPER_H
