//
// Created by mainf on 2024/11/30.
//

#ifndef SYMBOL_CALC_H
#define SYMBOL_CALC_H
#include <symbolic_analysis.h>
#include "base/matrix.h"

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus


void static_symbol_calc(SparseMatrix *A, INDEX_TYPE *l_count,INDEX_TYPE *u_count,
                     INDEX_TYPE *Rp_l, INDEX_TYPE *Ri_l,
                     INDEX_TYPE *Rp_u, INDEX_TYPE *Ri_u);

#ifdef __cplusplus
}
#endif //__cplusplus
#endif //SYMBOL_CALC_H
