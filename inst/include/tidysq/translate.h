#pragma once
#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err58-cpp"

#include <map>
#include "tidysq/ops/internal/util.h"

namespace tidysq {
    typedef std::unordered_map<LetterValue, std::unordered_map<LetterValue, std::unordered_map<LetterValue, LetterValue>>> CodonTable;

    CodonTable codon_table_1 = {
            {0u, {
                         {0u, {
                                      {0u, 8u}, {1u, 11u}, {2u, 8u}, {3u, 11u}
                         }},
                         {1u, {
                                      {0u, 16u}, {1u, 16u}, {2u, 16u}, {3u, 16u}
                         }},
                         {2u, {
                                      {0u, 14u}, {1u, 15u}, {2u, 14u}, {3u, 15u}
                         }},
                         {3u, {
                                      {0u, 7u}, {1u, 7u}, {2u, 10u}, {3u, 7u}
                         }}
            }},
            {1u, {
                         {0u, {
                                      {0u, 13u}, {1u, 6u}, {2u, 13u}, {3u, 6u}
                         }},
                         {1u, {
                                      {0u, 12u}, {1u, 12u}, {2u, 12u}, {3u, 12u}
                         }},
                         {2u, {
                                      {0u, 14u}, {1u, 14u}, {2u, 14u}, {3u, 14u}
                         }},
                         {3u, {
                                      {0u, 9u}, {1u, 9u}, {2u, 9u}, {3u, 9u}
                         }}
            }},
            {2u, {
                         {0u, {
                                      {0u, 3u}, {1u, 2u}, {2u, 3u}, {3u, 2u}
                         }},
                         {1u, {
                                      {0u, 0u}, {1u, 0u}, {2u, 0u}, {3u, 0u}
                         }},
                         {2u, {
                                      {0u, 5u}, {1u, 5u}, {2u, 5u}, {3u, 5u}
                         }},
                         {3u, {
                                      {0u, 17u}, {1u, 17u}, {2u, 17u}, {3u, 17u}
                         }}
            }},
            {3u, {
                         {0u, {
                                      {0u, 21u}, {1u, 19u}, {2u, 21u}, {3u, 19u}
                         }},
                         {1u, {
                                      {0u, 15u}, {1u, 15u}, {2u, 15u}, {3u, 15u}
                         }},
                         {2u, {
                                      {0u, 21u}, {1u, 1u}, {2u, 18u}, {3u, 1u}
                         }},
                         {3u, {
                                      {0u, 4u}, {1u, 9u}, {2u, 4u}, {3u, 9u}
                         }}
            }}
    };

    std::unordered_map<int, CodonTable> codon_diff_tables = {
            {2, {
                        {0u, {
                                     {2u, {
                                                  {0u, 21u}, {2u, 21u}
                                     }},
                                     {3u, {
                                                  {0u, 10u}
                                     }}
                        }},
                        {3u, {
                                     {2u, {
                                                  {0u, 18u}
                                     }}
                        }}
            }},
            {3, {
                        {0u, {
                                     {3u, {
                                                  {0u, 10u}
                                     }}
                        }},
                        {1u, {
                                     {2u, {
                                                  {0u, 31u}, {1u, 31u}
                                     }},
                                     {3u, {
                                                  {0u, 16u}, {1u, 16u}, {2u, 16u}, {3u, 16u}
                                     }}
                        }},
                        {3u, {
                                     {2u, {
                                                  {0u, 18u}
                                     }}
                        }}
            }},
            {4, {
                        {3u, {
                                     {2u, {
                                                  {0u, 18u}
                                     }}
                        }}
            }},
            {5, {
                        {0u, {
                                     {2u, {
                                                  {0u, 15u}, {2u, 15u}
                                     }},
                                     {3u, {
                                                  {0u, 10u}
                                     }}
                        }},
                        {3u, {
                                     {2u, {
                                                  {0u, 18u}
                                     }}
                        }}
            }},
            {6, {
                        {3u, {
                                     {0u, {
                                                  {0u, 13u}, {2u, 13u}
                                     }}
                        }}
            }},
            {9, {
                        {0u, {
                                     {0u, {
                                                  {0u, 11u}
                                     }},
                                     {2u, {
                                                  {0u, 15u}, {2u, 15u}
                                     }}
                        }},
                        {3u, {
                                     {2u, {
                                                  {0u, 18u}
                                     }}
                        }}
            }},
            {10, {
                         {3u, {
                                      {2u, {
                                                   {0u, 1u}
                                      }}
                         }}
            }},
            {12, {
                         {1u, {
                                      {3u, {
                                                   {2u, 15u}
                                      }}
                         }}
            }},
            {13, {
                         {0u, {
                                      {2u, {
                                                   {0u, 5u}, {2u, 5u}
                                      }},
                                      {3u, {
                                                   {0u, 10u}
                                      }}
                         }},
                         {3u, {
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {14, {
                         {0u, {
                                      {0u, {
                                                   {0u, 11u}
                                      }},
                                      {2u, {
                                                   {0u, 15u}, {2u, 15u}
                                      }}
                         }},
                         {3u, {
                                      {0u, {
                                                   {0u, 19u}
                                      }},
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {15, {
                         {3u, {
                                      {0u, {
                                                   {2u, 13u}
                                           }}
                              }}
            }},
            {16, {
                         {3u, {
                                      {0u, {
                                                   {2u, 9u}
                                      }}
                         }}
            }},
            {21, {
                         {0u, {
                                      {0u, {
                                                   {0u, 11u}
                                      }},
                                      {2u, {
                                                   {0u, 15u}, {2u, 15u}
                                      }},
                                      {3u, {
                                                   {0u, 10u}
                                      }}
                         }},
                         {3u, {
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {22, {
                         {3u, {
                                      {0u, {
                                                   {2u, 9u}
                                      }},
                                      {1u, {
                                                   {0u, 21u}
                                      }}
                         }}
            }},
            {23, {
                         {3u, {
                                      {3u, {
                                                   {0u, 21u}
                                      }}
                         }}
            }},
            {24, {
                         {0u, {
                                      {2u, {
                                                   {0u, 15u}, {2u, 8u}
                                      }}
                         }},
                         {3u, {
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {25, {
                         {3u, {
                                      {2u, {
                                                   {0u, 5u}
                                      }}
                         }}
            }},
            {26, {
                         {1u, {
                                      {3u, {
                                                   {2u, 0u}
                                      }}
                         }}
            }},
            {27, {
                         {3u, {
                                      {0u, {
                                                   {0u, 13u}, {2u, 13u}
                                      }},
                         }}
            }},
            {29, {
                         {3u, {
                                      {0u, {
                                                   {0u, 19u}, {2u, 19u}
                                      }}
                         }}
            }},
            {30, {
                         {3u, {
                                      {0u, {
                                                   {0u, 3u}, {2u, 3u}
                                      }}
                         }}
            }},
            {31, {
                         {3u, {
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {33, {
                         {0u, {
                                      {2u, {
                                                   {0u, 15u}, {2u, 8u}
                                      }}
                         }},
                         {3u, {
                                      {0u, {
                                                   {0u, 19u}
                                      }},
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }}
    };

    std::unordered_map<int, CodonTable> amb_codon_diff_tables = {
            {27, {
                         {3u, {
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {28, {
                         {3u, {
                                      {0u, {
                                                   {0u, 13u}, {2u, 13u}
                                      }},
                                      {2u, {
                                                   {0u, 18u}
                                      }}
                         }}
            }},
            {31, {
                         {3u, {
                                      {0u, {
                                                   {0u, 3u}, {2u, 3u}
                                      }}
                         }}
            }}
    };

    LetterValue codon_table(int table,
                            const LetterValue &codon_1,
                            const LetterValue &codon_2,
                            const LetterValue &codon_3,
                            const bool &interpret_as_stop) {
        // Some tables are actually identical to some others
        if (table == 7) table = 4;
        if (table == 8) table = 1;
        if (table == 11) table = 1;
        // If this is non-standard table, then we have this kind of @Override
        if (table != 1) {
            // The only way to include ambiguous translation tables (27, 28, 31)
            // I don't really like this solution, maybe we should just drop the support for that
            if (table != 27 && table != 28 && table != 31 && !interpret_as_stop) {
                auto amb_codon_diff_table = amb_codon_diff_tables[table];
                if (amb_codon_diff_table.count(codon_1) > 0 &&
                        amb_codon_diff_table[codon_1].count(codon_2) > 0 &&
                        amb_codon_diff_table[codon_1][codon_2].count(codon_3) > 0) {
                    auto amino_acid = amb_codon_diff_table[codon_1][codon_2][codon_3];
                    return amino_acid;
                }
            }
            // Then find correct table of differences
            auto codon_diff_table = codon_diff_tables[table];
            if (codon_diff_table.count(codon_1) > 0 &&
                    codon_diff_table[codon_1].count(codon_2) > 0 &&
                    codon_diff_table[codon_1][codon_2].count(codon_3) > 0) {
                auto amino_acid = codon_diff_table[codon_1][codon_2][codon_3];
                // TODO: handle case if (amino_acid == 31u) (i.e. NA_letter)
                return amino_acid;
            }
        }
        return codon_table_1[codon_1][codon_2][codon_3];
    }

    template<InternalType INTERNAL>
    Sequence<INTERNAL> translate(const Sequence<INTERNAL> &sequence,
                                 const int &table,
                                 const bool &interpret_as_stop,
                                 const AlphSize &input_alph_size,
                                 const AlphSize &output_alph_size) {
        LenSq sequence_length = sequence.originalLength() / 3;
        Sequence<INTERNAL> ret = internal::reserve_space_for_packed<INTERNAL>(sequence_length, output_alph_size);

        if (sequence_length > 0) {
            auto input_it = sequence.cbegin(input_alph_size);
            auto output_it = ret.begin(output_alph_size);
            while (input_it < sequence.cend(input_alph_size) - 2) {
                auto codon_1 = *input_it++;
                auto codon_2 = *input_it++;
                auto codon_3 = *input_it++;
                output_it.assign(codon_table(table, codon_1, codon_2, codon_3, interpret_as_stop));
                ++output_it;
            }
        }
        return ret;
    }

    template<InternalType INTERNAL>
    Sq<INTERNAL> translate(const Sq<INTERNAL> &sq,
                           const int &table,
                           const bool &interpret_as_stop) {
        const Alphabet& alph = sq.alphabet();
        Sq<INTERNAL> ret(sq.length(), Alphabet(AMI_BSC, alph.NA_letter()));

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = translate(sq[i].get(), table, interpret_as_stop,
                    alph.alphabet_size(), ret.alphabet().alphabet_size());
        }
        return ret;
    }
}

#pragma clang diagnostic pop