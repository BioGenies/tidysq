#pragma once 

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/unpack.h"
#include "tidysq/ops/pack.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename PROTO_TMP>
        class OperationApplyRFunction : public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                                                                        Rcpp::List, SEXP> {
            const Alphabet &alphabet_;
            const Rcpp::Function &function_;

        public:
            explicit OperationApplyRFunction(const Alphabet &alphabet,
                                             const Rcpp::Function &function) :
                    alphabet_(alphabet),
                    function_(function) {};

            inline Rcpp::List initialize_vector_out(const Sq<INTERNAL_IN> &sq, const LenSq from, const LenSq to) override {
                auto ret = Rcpp::List::create();
                for (LenSq i = 0; i < to - from; i++) {
                    ret.push_back(Rcpp::IntegerVector{});
                }
                return ret;
            }

            inline SEXP initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return {};
            }

            inline void operator() (const Sequence<INTERNAL_IN> &sequence,
                                    SEXP &ret_elem) override {
                auto unpacked = tidysq::unpack<INTERNAL_IN, RCPP_IT, PROTO_TMP>(sequence, alphabet_);
                ret_elem = function_(unpacked.content());
            }

            inline SEXP operator() (const Sequence<INTERNAL_IN> &sequence) override {
                SEXP ret_elem = initialize_element_out(sequence);
                operator()(sequence, ret_elem);
                return ret_elem;
            }
        };
    }

    template<typename INTERNAL_IN, typename PROTO_TMP>
    Rcpp::List apply_R_function(const Sq<INTERNAL_IN> &sq,
                                const Rcpp::Function &function) {
        return sqapply(sq, ops::OperationApplyRFunction<INTERNAL_IN, PROTO_TMP>(sq.alphabet(), function));
    }
}