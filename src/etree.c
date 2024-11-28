//
// Created by mainf on 2024/11/28.
//

#include "etree.h"

/*
 *  Implementation of disjoint set union routines.
 *  Elements are integers in 0..n-1, and the
 *  names of the sets themselves are of type int.
 *
 *  Calls are:
 *  initialize_disjoint_sets (n) initial call.
 *  s = make_set (i)             returns a set containing only i.
 *  s = link (t, u)		 returns s = t union u, destroying t and u.
 *  s = find (i)		 return name of set containing i.
 *  finalize_disjoint_sets 	 final call.
 *
 *  This implementation uses path compression but not weighted union.
 *  See Tarjan's book for details.
 *  John Gilbert, CMI, 1987.
 *
 *  Implemented path-halving by XSL 07/05/95.
 */

static
INDEX_TYPE *mxCallocInt(INDEX_TYPE n) {
    register INDEX_TYPE i;
    INDEX_TYPE *buf;

    buf = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    if (!buf) {
        LOG_ERROR("lu_malloc fails for buf in mxCallocInt()");
    }
    for (i = 0; i < n; i++) buf[i] = 0;
    return (buf);
}

static
void initialize_disjoint_sets(
    INDEX_TYPE n,
    INDEX_TYPE **pp
) {
    (*pp) = mxCallocInt(n);
}


static
INDEX_TYPE make_set(
    INDEX_TYPE i,
    INDEX_TYPE *pp
) {
    pp[i] = i;
    return i;
}


static
INDEX_TYPE link(
    INDEX_TYPE s,
    INDEX_TYPE t,
    INDEX_TYPE *pp
) {
    pp[s] = t;
    return t;
}


/* PATH HALVING */
static INDEX_TYPE find(INDEX_TYPE i, INDEX_TYPE *pp) {
    register INDEX_TYPE p, gp;

    p = pp[i];
    gp = pp[p];
    while (gp != p) {
        pp[i] = gp;
        i = gp;
        p = pp[i];
        gp = pp[p];
    }
    return (p);
}

#if 0
/* PATH COMPRESSION */
static
int find (
	int i
	)
{
	if (pp[i] != i)
		pp[i] = find (pp[i]);
	return pp[i];
}
#endif

static
void finalize_disjoint_sets(
    INDEX_TYPE *pp
) {
    lu_free(pp);
}

/*
 * Symmetric elimination tree
 */
INDEX_TYPE
sp_symetree(
    INDEX_TYPE *acolst, INDEX_TYPE *acolend, /* column starts and ends past 1 */
    INDEX_TYPE *arow, /* row indices of A */
    INDEX_TYPE n, /* dimension of A */
    INDEX_TYPE *parent /* parent in elim tree */
) {
    INDEX_TYPE *root; /* root of subtree of etree 	*/
    INDEX_TYPE rset, cset;
    INDEX_TYPE row, col;
    INDEX_TYPE rroot;
    INDEX_TYPE p;
    INDEX_TYPE *pp;

    root = mxCallocInt(n);
    initialize_disjoint_sets(n, &pp);

    for (col = 0; col < n; col++) {
        cset = make_set(col, pp);
        root[cset] = col;
        parent[col] = -1; /* Matlab */
        for (p = acolst[col]; p < acolend[col]; p++) {
            row = arow[p];
            if (row >= col) continue;
            rset = find(row, pp);
            rroot = root[rset];
            if (rroot != col) {
                parent[rroot] = col;
                cset = link(cset, rset, pp);
                root[cset] = col;
            }
        }
    }
    lu_free(root);
    finalize_disjoint_sets(pp);
    return 0;
} /* SP_SYMETREE */
