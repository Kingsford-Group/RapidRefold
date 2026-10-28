#ifndef RAPIDREFOLD_MFE_INCREMENTAL_H
#define RAPIDREFOLD_MFE_INCREMENTAL_H

/**
 *  @file     mfe_incremental.h
 *  @brief    Incremental MFE recomputation after point mutations
 *
 *  After an initial MFE computation via vrna_mfe(), a point mutation
 *  only invalidates DP subproblems whose interval [i,j] contains the
 *  mutated position. This module re-fills only those entries, avoiding
 *  a full O(n^3) recomputation.
 *
 *  Methods are composable via orthogonal flags:
 *
 *  Setup axis (how to update ptype/hc after mutation):
 *    0                    -- O(n^2) full ptype rebuild + vrna_hc_init
 *    VRNA_MFE_INC_LAZY   -- O(n) update only entries involving pos
 *
 *  DP fill axis (which fML rows to recompute):
 *    0                    -- recompute all rows i <= pos+1
 *    VRNA_MFE_INC_SPARSE -- skip rows where no c/fML value changed
 *
 *  Bifurcation axis (how to compute the fML split decomposition):
 *    0                    -- O(n) full scan per entry
 *    VRNA_MFE_INC_TOPK   -- O(K) lookup from pre-cached split points
 */

#include <ViennaRNA/fold_compound.h>

/**
 *  @brief  Cache of top-K split points and interior loop candidates
 *
 *  Stores the K best candidates for each (i,j) entry from the
 *  initial fold. Used by VRNA_MFE_INC_TOPK to replace O(n) or
 *  O(MAXLOOP^2) scans with O(K) lookups.
 */
typedef struct {
  unsigned short  *fml_splits;    /**< top-K splits for fML bifurcation */
  unsigned short  *mb_splits;     /**< top-K splits for MB loop closing pair */
  unsigned short  *int_k;         /**< top-K interior loop inner pair 5' pos */
  unsigned short  *int_l;         /**< top-K interior loop inner pair 3' pos */
  unsigned char   *fml_count;     /**< valid splits stored for fML at (i,j) */
  unsigned char   *mb_count;      /**< valid splits stored for MB at (i,j) */
  unsigned char   *int_count;     /**< valid interior loop candidates at (i,j) */
  unsigned int    K;              /**< max candidates per entry */
  unsigned int    n;              /**< sequence length */
  int             *indx;          /**< pointer to fc->jindx */
} vrna_topk_cache_t;


/**
 *  @brief  Build a top-K cache from a folded compound
 *
 *  The fold compound must already have filled DP matrices from
 *  a prior call to vrna_mfe(). Scans the fML and c matrices to
 *  record the K best bifurcation splits and interior loop pairs.
 *
 *  @param  fc  Fold compound with existing DP matrices
 *  @param  K   Number of top candidates to cache per entry
 *  @return     Allocated cache, or NULL on error. Free with vrna_topk_cache_free().
 */
vrna_topk_cache_t *
vrna_topk_cache_from_fold(vrna_fold_compound_t *fc,
                          unsigned int         K);


/**
 *  @brief  Free a top-K cache
 */
void
vrna_topk_cache_free(vrna_topk_cache_t *cache);


/* ================================================================
 * Composable flags
 * ================================================================ */

/** @brief Setup: O(n) ptype/hc update instead of O(n^2) full rebuild */
#define VRNA_MFE_INC_LAZY     0x01u

/** @brief DP fill: skip fML rows where no c/fML value changed */
#define VRNA_MFE_INC_SPARSE   0x02u

/** @brief Bifurcation: use top-K cached split points instead of O(n) scan */
#define VRNA_MFE_INC_TOPK     0x04u

/** @brief All exact optimizations (lazy setup + sparse fill) */
#define VRNA_MFE_INC_ALL      (VRNA_MFE_INC_LAZY | VRNA_MFE_INC_SPARSE)


/**
 *  @brief  Recompute MFE after a single-nucleotide mutation (composable)
 *
 *  The fold compound must already have filled DP matrices from
 *  a prior call to vrna_mfe(). Mutates position pos (1-based) to
 *  new_base, updates internal data, then incrementally recomputes
 *  only the affected DP entries.
 *
 *  Exact combinations (produce identical results to full recomputation):
 *    0, LAZY, SPARSE, LAZY|SPARSE
 *
 *  Heuristic combinations (may differ from full recomputation):
 *    Any combination including TOPK
 *
 *  @param  fc            Fold compound with existing DP matrices
 *  @param  pos           1-based mutation position
 *  @param  new_base      The new nucleotide character
 *  @param  flags         Composable flags (OR of VRNA_MFE_INC_*)
 *  @param  structure     Output buffer for the new structure (may be NULL)
 *  @param  cache         Top-K cache (required if TOPK flag set, else NULL)
 *  @return               New MFE in kcal/mol, or INF/100 on error
 */
float
vrna_mfe_incremental_ex(vrna_fold_compound_t *fc,
                        unsigned int         pos,
                        char                 new_base,
                        unsigned int         flags,
                        char                 *structure,
                        vrna_topk_cache_t    *cache);


/**
 *  @brief  Recompute MFE after a single-nucleotide mutation (base method)
 *
 *  Convenience wrapper: calls vrna_mfe_incremental_ex() with flags=0.
 */
float
vrna_mfe_incremental(vrna_fold_compound_t *fc,
                     unsigned int         pos,
                     char                 new_base,
                     char                 *structure);

#endif
