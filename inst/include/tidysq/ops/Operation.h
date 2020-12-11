#pragma once

#include "tidysq/Alphabet.h"

namespace tidysq::ops {
    template<typename VECTOR_IN, typename ELEMENT_IN,
             typename VECTOR_OUT, typename ELEMENT_OUT>
    class OperationVectorToVector {
    public:
        virtual bool may_return_early(const VECTOR_IN &vector_in) {
            return false;
        }

        virtual VECTOR_OUT return_early(const VECTOR_IN &vector_in) {
            return  initialize_vector_out(vector_in);
        }

        virtual VECTOR_OUT initialize_vector_out(const VECTOR_IN &vector_in, const LenSq from, const LenSq to) = 0;

        inline virtual VECTOR_OUT initialize_vector_out(const VECTOR_IN &vector_in) {
            return initialize_vector_out(vector_in, 0, vector_in.size());
        }

        virtual ELEMENT_OUT initialize_element_out(const ELEMENT_IN &element_in) = 0;

        virtual void operator() (const ELEMENT_IN &element_in, ELEMENT_OUT &element_out) = 0;

        //TODO: issue #57
        virtual ELEMENT_OUT operator() (const ELEMENT_IN &element_in) {
            ELEMENT_OUT element_out = initialize_element_out(element_in);
            operator()(element_in, element_out);
            return element_out;
        }
    };

    template<typename INTERNAL_IN, typename INTERNAL_OUT>
    class OperationSqToSq :
            public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                                           Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>> {
    public:
        inline Sq<INTERNAL_OUT> initialize_vector_out(const Sq<INTERNAL_IN> &sq_in, LenSq from, LenSq to) override {
            return Sq<INTERNAL_OUT>(to - from, map_alphabet(sq_in.alphabet()));
        }

        inline Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence_in) override {
            return Sequence<INTERNAL_OUT>(sequence_in.size(), sequence_in.original_length());
        }

        [[nodiscard]] virtual inline Alphabet map_alphabet(const Alphabet &alphabet_in) const {
            return alphabet_in;
        }
    };
}
