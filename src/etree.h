//
// Created by mainf on 2024/11/28.
//

#ifndef ETREE_H
#define ETREE_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <base/malloc.h>
#include <base/matrix.h>

/**
 * Symmetric elimination tree
 **/
INDEX_TYPE sp_symetree(
    INDEX_TYPE *acolst,
    INDEX_TYPE *acolend, /* column starts and ends past 1 */
    INDEX_TYPE *arow, /* row indices of A */
    INDEX_TYPE n, /* dimension of A */
    INDEX_TYPE *parent /* parent in elim tree */
);

#ifdef __cplusplus
}
#endif

#endif //ETREE_H
