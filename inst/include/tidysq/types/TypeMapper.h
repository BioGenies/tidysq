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
    /*  InternalTypeMapper<INTERNAL> keeps definition of the following types:
     *      - SequenceElementType := element of Sequence<INTERNAL>::content_
     *      - SequenceContentType := decltype(Sequence<INTERNAL>::content_)
     *      - SequenceAccessType := decltype(typename SequenceContentType::operator[](LenSq))
     *      - SequenceConstAccessType := decltype(typename SequenceContentType::operator[](LenSq) const)
     *
     *      - SqContentType := decltype(Sq<INTERNAL>::content_)
     *      - SqElementRealType := element of Sq<INTERNAL>::content_
     *      - SqElementDestType := Sequence<INTERNAL>
     *      - SqAccessType := decltype(typename SqContentType::operator[](LenSq))
     *      - SqConstAccessType := decltype(typename SqContentType::operator[](LenSq) const)
     */

    template<>
    struct InternalTypeMapper<STD> {
        typedef ElemPacked                              SequenceElementType;
        typedef std::vector<ElemPacked>                 SequenceContentType;
        typedef         SequenceElementType &           SequenceAccessType;
        typedef const   SequenceElementType &           SequenceConstAccessType;

        typedef std::vector<Sequence<STD>>              SqContentType;
        typedef Sequence<STD>                           SqElementType;
        typedef         Sequence<STD> &                 SqAccessType;
        typedef const   Sequence<STD> &                 SqConstAccessType;
    };

    template<>
    struct InternalTypeMapper<RCPP> {
        typedef ElemPacked                              SequenceElementType;
        typedef Rcpp::RawVector                         SequenceContentType;
        typedef Rcpp::RawVector::Proxy                  SequenceAccessType;
        typedef Rcpp::RawVector::const_Proxy            SequenceConstAccessType;

        typedef Rcpp::List                              SqContentType;
        typedef Sequence<RCPP>                          SqElementType;
        typedef         Rcpp::List::Proxy               SqAccessType;
        typedef const   Rcpp::List::const_Proxy         SqConstAccessType;
    };

    template<ProtoType PROTO>
    struct ProtoTypeMapper;

    template<>
    struct ProtoTypeMapper<RAWS> {
        typedef Rcpp::RawVector                         ProtoSequenceRcppContentType;
    };

    template<InternalType INTERNAL, ProtoType PROTO>
    struct TypeMapper;
    
    template<>
    struct TypeMapper<STD, RAWS> {
        typedef ElemRaws                                ProtoSequenceElementType;
        typedef std::vector<ProtoSequenceElementType>   ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &      ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &      ProtoSequenceConstAccessType;

        typedef std::vector<ProtoSequence<STD, RAWS>>   ProtoSqContentType;
        typedef ProtoSequence<STD, RAWS>                ProtoSqElementType;
        typedef         ProtoSequence<STD, RAWS> &      ProtoSqAccessType;
        typedef const   ProtoSequence<STD, RAWS> &      ProtoSqConstAccessType;

    };

    template<>
    struct TypeMapper<RCPP, RAWS> {
        typedef ElemRaws                                ProtoSequenceElementType;
        typedef Rcpp::RawVector                         ProtoSequenceContentType;
        typedef Rcpp::RawVector::Proxy                  ProtoSequenceAccessType;
        typedef Rcpp::RawVector::const_Proxy            ProtoSequenceConstAccessType;

        typedef Rcpp::List                              ProtoSqContentType;
        typedef ProtoSequence<RCPP, RAWS>               ProtoSqElementType;
        typedef         Rcpp::List::Proxy               ProtoSqAccessType;
        typedef const   Rcpp::List::const_Proxy         ProtoSqConstAccessType;
    };

    template<>
    struct TypeMapper<STD, INTS> {
        typedef ElemInts                                ProtoSequenceElementType;
        typedef std::vector<ProtoSequenceElementType>   ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &      ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &      ProtoSequenceConstAccessType;

        typedef std::vector<ProtoSequence<STD, INTS>>   ProtoSqContentType;
        typedef ProtoSequence<STD, INTS>                ProtoSqElementType;
        typedef         ProtoSequence<STD, INTS> &      ProtoSqAccessType;
        typedef const   ProtoSequence<STD, INTS> &      ProtoSqConstAccessType;

    };

    template<>
    struct TypeMapper<RCPP, INTS> {
        typedef ElemInts                                ProtoSequenceElementType;
        typedef Rcpp::IntegerVector                     ProtoSequenceContentType;
        typedef Rcpp::IntegerVector::Proxy              ProtoSequenceAccessType;
        typedef Rcpp::IntegerVector::const_Proxy        ProtoSequenceConstAccessType;

        typedef Rcpp::List                              ProtoSqContentType;
        typedef ProtoSequence<RCPP, INTS>               ProtoSqElementType;
        typedef         Rcpp::List::Proxy               ProtoSqAccessType;
        typedef const   Rcpp::List::const_Proxy         ProtoSqConstAccessType;
    };

    template<>
    struct TypeMapper<STD, STRINGS> {
        typedef ElemStrings                             ProtoSequenceElementType;
        typedef std::vector<ProtoSequenceElementType>   ProtoSequenceContentType;
        typedef         ProtoSequenceElementType &      ProtoSequenceAccessType;
        typedef const   ProtoSequenceElementType &      ProtoSequenceConstAccessType;

        typedef std::vector<ProtoSequence<STD, RAWS>>   ProtoSqContentType;
        typedef ProtoSequence<STD, RAWS>                ProtoSqElementType;
        typedef         ProtoSequence<STD, RAWS> &      ProtoSqAccessType;
        typedef const   ProtoSequence<STD, RAWS> &      ProtoSqConstAccessType;

    };

    template<>
    struct TypeMapper<RCPP, STRINGS> {
        typedef ElemStrings                             ProtoSequenceElementType;
        typedef Rcpp::StringVector                      ProtoSequenceContentType;
        typedef Rcpp::StringVector::Proxy               ProtoSequenceAccessType;
        typedef Rcpp::StringVector::const_Proxy         ProtoSequenceConstAccessType;

        typedef Rcpp::List                              ProtoSqContentType;
        typedef ProtoSequence<RCPP, STRINGS>            ProtoSqElementType;
        typedef         Rcpp::List::Proxy               ProtoSqAccessType;
        typedef const   Rcpp::List::const_Proxy         ProtoSqConstAccessType;
    };

}

#endif //TIDYSQ_TYPEMAPPER_H
