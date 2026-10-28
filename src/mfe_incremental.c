/*
 *  Incremental MFE recomputation after point mutations.
 *
 *  When a single nucleotide changes at position k, only DP entries
 *  c[i][j], fML[i][j], fM1[i][j] where i <= k <= j are affected.
 *  This module updates the sequence data, then re-runs the DP fill
 *  selectively for those entries, followed by a full f5 recomputation.
 *
 *  The fML computation uses the same E_ml_stems_fast path as the
 *  original algorithm (via vrna_E_ml_stems_fast) with proper rotating
 *  auxiliary arrays, ensuring exact consistency with the original DP.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/structures/dotbracket.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/mfe/exterior.h"
#include "ViennaRNA/mfe/internal.h"
#include "ViennaRNA/mfe/multibranch.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/backtrack/global.h"

#include "mfe_incremental.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define ADD_OR_INF(a, b) (((a) != INF) && ((b) != INF) ? (a) + (b) : INF)

extern int
vrna_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli);


static INLINE int
incremental_decompose_pair(vrna_fold_compound_t *fc,
                           int                  i,
                           int                  j);

static INLINE int
incremental_decompose_pair_topk(vrna_fold_compound_t *fc,
                                int                  i,
                                int                  j,
                                vrna_topk_cache_t    *cache);

static INLINE int
topk_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     int                  *Fmi,
                     int                  *DMLi,
                     vrna_topk_cache_t    *cache);


/*
 *  Update sequence encodings and pair-type table at a single position.
 */
static void
update_sequence_at_pos(vrna_fold_compound_t *fc,
                       unsigned int         pos,
                       char                 new_base)
{
  vrna_md_t     *md = &(fc->params->model_details);
  unsigned int  n   = fc->length;

  new_base = toupper(new_base);

  fc->sequence[pos - 1] = new_base;

  if (fc->nucleotides) {
    fc->nucleotides[0].string[pos - 1] = new_base;

    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);
    short alias_enc  = md->alias[simple_enc];

    fc->nucleotides[0].encoding[pos] = alias_enc;

    if (pos < fc->nucleotides[0].length)
      fc->nucleotides[0].encoding5[pos + 1] = alias_enc;

    if (pos > 1)
      fc->nucleotides[0].encoding3[pos - 1] = alias_enc;
  }

  {
    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);
    short alias_enc  = md->alias[simple_enc];

    fc->sequence_encoding[pos] = alias_enc;

    if (pos == 1)
      fc->sequence_encoding[n + 1] = alias_enc;
    if (pos == n)
      fc->sequence_encoding[0] = alias_enc;
  }

  {
    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);

    fc->sequence_encoding2[pos] = simple_enc;

    if (pos == 1)
      fc->sequence_encoding2[n + 1] = simple_enc;
  }

  if (fc->encoding5 && pos < n)
    fc->encoding5[pos + 1] = fc->sequence_encoding[pos];

  if (fc->encoding3 && pos > 1)
    fc->encoding3[pos - 1] = fc->sequence_encoding[pos];

  if (fc->ptype) {
    free(fc->ptype);
    fc->ptype = vrna_ptypes(fc->sequence_encoding2, md);
  }
}


/*
 *  Update only the ptype entries involving position pos.
 *  O(n) instead of O(n^2) for the full vrna_ptypes() rebuild.
 */
static void
update_ptype_at_pos(vrna_fold_compound_t *fc,
                    unsigned int         pos)
{
  unsigned int  i, j, n;
  int           *indx;
  short         *S2;
  vrna_md_t     *md;
  char          *ptype;

  n     = fc->length;
  indx  = fc->jindx;
  S2    = fc->sequence_encoding2;
  md    = &(fc->params->model_details);
  ptype = fc->ptype;

  if (!ptype)
    return;

  for (i = 1; i < pos; i++) {
    if ((pos - i) > (unsigned int)md->max_bp_span)
      continue;
    if ((pos - i - 1) < (unsigned int)md->min_loop_size)
      continue;

    int type = md->pair[S2[i]][S2[pos]];

    if (md->noLP) {
      int otype = 0, ntype = 0;
      if (i > 1 && pos < n)
        otype = md->pair[S2[i - 1]][S2[pos + 1]];
      if (i + 1 < pos && (pos - i - 2) >= (unsigned int)md->min_loop_size)
        ntype = md->pair[S2[i + 1]][S2[pos - 1]];

      if (!otype && !ntype)
        type = 0;
    }

    ptype[indx[pos] + i] = (char)type;
  }

  for (j = pos + 1; j <= n; j++) {
    if ((j - pos) > (unsigned int)md->max_bp_span)
      continue;
    if ((j - pos - 1) < (unsigned int)md->min_loop_size)
      continue;

    int type = md->pair[S2[pos]][S2[j]];

    if (md->noLP) {
      int otype = 0, ntype = 0;
      if (pos > 1 && j < n)
        otype = md->pair[S2[pos - 1]][S2[j + 1]];
      if (pos + 1 < j && (j - pos - 2) >= (unsigned int)md->min_loop_size)
        ntype = md->pair[S2[pos + 1]][S2[j - 1]];

      if (!otype && !ntype)
        type = 0;
    }

    ptype[indx[j] + pos] = (char)type;
  }

  if (md->noLP) {
    unsigned int neighbors[] = { pos > 1 ? pos - 1 : 0, pos < n ? pos + 1 : 0 };
    for (int nb = 0; nb < 2; nb++) {
      unsigned int p = neighbors[nb];
      if (p == 0) continue;

      for (j = p + 1; j <= n; j++) {
        if ((j - p) > (unsigned int)md->max_bp_span) continue;
        if ((j - p - 1) < (unsigned int)md->min_loop_size) continue;

        int type = md->pair[S2[p]][S2[j]];
        if (type) {
          int otype = 0, ntype = 0;
          if (p > 1 && j < n)
            otype = md->pair[S2[p - 1]][S2[j + 1]];
          if (p + 1 < j && (j - p - 2) >= (unsigned int)md->min_loop_size)
            ntype = md->pair[S2[p + 1]][S2[j - 1]];
          if (!otype && !ntype) type = 0;
        }
        ptype[indx[j] + p] = (char)type;
      }

      for (i = 1; i < p; i++) {
        if ((p - i) > (unsigned int)md->max_bp_span) continue;
        if ((p - i - 1) < (unsigned int)md->min_loop_size) continue;

        int type = md->pair[S2[i]][S2[p]];
        if (type) {
          int otype = 0, ntype = 0;
          if (i > 1 && p < n)
            otype = md->pair[S2[i - 1]][S2[p + 1]];
          if (i + 1 < p && (p - i - 2) >= (unsigned int)md->min_loop_size)
            ntype = md->pair[S2[i + 1]][S2[p - 1]];
          if (!otype && !ntype) type = 0;
        }
        ptype[indx[p] + i] = (char)type;
      }
    }
  }
}


/*
 *  Update hc->mx entries involving position pos.  O(n).
 */
static void
update_hc_at_pos(vrna_fold_compound_t *fc,
                 unsigned int         pos)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;
  short         *S2;
  vrna_md_t     *md;

  n   = fc->length;
  hc  = fc->hc;
  S2  = fc->sequence_encoding2;
  md  = &(fc->params->model_details);

  if (!hc || !hc->mx)
    return;

  unsigned int affected[] = { pos > 1 ? pos - 1 : 0, pos, pos < n ? pos + 1 : 0 };
  for (int a = 0; a < 3; a++) {
    unsigned int p = affected[a];
    if (p == 0) continue;

    for (i = 1; i <= n; i++) {
      if (i == p) continue;
      unsigned int lo = (i < p) ? i : p;
      unsigned int hi = (i < p) ? p : i;

      unsigned char constraint = VRNA_CONSTRAINT_CONTEXT_NONE;

      if (((hi - lo) < (unsigned int)md->max_bp_span) &&
          ((fc->strand_number[lo] != fc->strand_number[hi]) ||
           ((hi - lo) > (unsigned int)md->min_loop_size))) {
        int type = md->pair[S2[lo]][S2[hi]];
        switch (type) {
          case 0:
            break;
          case 3:
          case 4:
            if (md->noGU)
              break;
            else if (md->noGUclosure) {
              constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              constraint &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP |
                              VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
              break;
            }
            /* fallthrough */
          default:
            constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            break;
        }

        if (md->noLP && constraint) {
          unsigned char can_stack = VRNA_CONSTRAINT_CONTEXT_NONE;
          if (lo > 1 && hi < n &&
              ((hi - lo + 2) < (unsigned int)md->max_bp_span ||
               fc->strand_number[lo - 1] != fc->strand_number[hi + 1]) &&
              md->pair[S2[lo - 1]][S2[hi + 1]])
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          if (lo + 2 < hi &&
              ((hi - lo - 2) > (unsigned int)md->min_loop_size ||
               fc->strand_number[lo + 1] != fc->strand_number[hi - 1]) &&
              md->pair[S2[lo + 1]][S2[hi - 1]])
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          constraint &= can_stack;
        }
      }

      hc->mx[n * lo + hi] = constraint;
      hc->mx[n * hi + lo] = constraint;
    }
  }

  for (int a = 0; a < 3; a++) {
    unsigned int p = affected[a];
    if (p == 0) continue;
    hc->mx[n * p + p] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }

  for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--)
    hc->up_ext[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ?
                    1 + hc->up_ext[i + 1] : 0;

  for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)
    hc->up_hp[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ?
                   1 + hc->up_hp[i + 1] : 0;

  for (hc->up_int[n + 1] = 0, i = n; i > 0; i--)
    hc->up_int[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ?
                    1 + hc->up_int[i + 1] : 0;

  for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)
    hc->up_ml[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ?
                   1 + hc->up_ml[i + 1] : 0;
}


/*
 *  Lazy sequence update: encodings only, no ptype rebuild.
 */
static void
update_sequence_at_pos_lazy(vrna_fold_compound_t *fc,
                            unsigned int         pos,
                            char                 new_base)
{
  vrna_md_t     *md = &(fc->params->model_details);
  unsigned int  n   = fc->length;

  new_base = toupper(new_base);

  fc->sequence[pos - 1] = new_base;

  if (fc->nucleotides) {
    fc->nucleotides[0].string[pos - 1] = new_base;
    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);
    short alias_enc  = md->alias[simple_enc];
    fc->nucleotides[0].encoding[pos] = alias_enc;
    if (pos < fc->nucleotides[0].length)
      fc->nucleotides[0].encoding5[pos + 1] = alias_enc;
    if (pos > 1)
      fc->nucleotides[0].encoding3[pos - 1] = alias_enc;
  }

  {
    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);
    short alias_enc  = md->alias[simple_enc];
    fc->sequence_encoding[pos] = alias_enc;
    if (pos == 1) fc->sequence_encoding[n + 1] = alias_enc;
    if (pos == n) fc->sequence_encoding[0] = alias_enc;
  }

  {
    short simple_enc = (short)vrna_nucleotide_encode(new_base, md);
    fc->sequence_encoding2[pos] = simple_enc;
    if (pos == 1) fc->sequence_encoding2[n + 1] = simple_enc;
  }

  if (fc->encoding5 && pos < n)
    fc->encoding5[pos + 1] = fc->sequence_encoding[pos];

  if (fc->encoding3 && pos > 1)
    fc->encoding3[pos - 1] = fc->sequence_encoding[pos];
}


static INLINE int
mb_loop_from_fml(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j)
{
  short         *S2;
  unsigned int  tt, k, n;
  int           decomp, e, *fML, *indx, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n             = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  S2            = fc->sequence_encoding2;
  indx          = fc->jindx;
  fML           = fc->matrices->fML;
  hc            = fc->hc;

  decomp = INF;

  if (!(hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP))
    return INF;

  if (hc->eval_mb &&
      !hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc))
    return INF;

  for (k = i + 2; k <= j - 2; k++) {
    int fml_left  = fML[indx[k] + i + 1];
    int fml_right = fML[indx[j - 1] + k + 1];
    e = ADD_OR_INF(fml_left, fml_right);
    decomp = MIN2(decomp, e);
  }

  if (decomp != INF) {
    tt = vrna_get_ptype_md(S2[j], S2[i], md);

    if (md->noGUclosure && ((tt == 3) || (tt == 4)))
      return INF;

    switch (dangle_model) {
      case 2:
        decomp += vrna_E_multibranch_stem(tt, S2[j - 1], S2[i + 1], P) +
                  P->MLclosing;
        break;
      default:
        decomp += vrna_E_multibranch_stem(tt, -1, -1, P) +
                  P->MLclosing;
        break;
    }
  }

  return decomp;
}


static INLINE int
incremental_decompose_pair(vrna_fold_compound_t *fc,
                           int                  i,
                           int                  j)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, dangle_model;

  n             = fc->length;
  dangle_model  = fc->params->model_details.dangles;
  hc_decompose  = fc->hc->mx[n * i + j];
  e             = INF;

  if (hc_decompose) {
    new_c = INF;

    energy  = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
    new_c   = MIN2(new_c, energy);

    energy  = mb_loop_from_fml(fc, i, j);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      energy  = vrna_mfe_multibranch_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    energy  = vrna_mfe_internal(fc, i, j);
    new_c   = MIN2(new_c, energy);

    e = new_c;

    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->c); c++) {
        if (fc->aux_grammar->c[c].cb) {
          energy  = fc->aux_grammar->c[c].cb(fc, i, j, fc->aux_grammar->c[c].data);
          e       = MIN2(e, energy);
        }
      }
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (e != INF))
      e -= fc->pscore[fc->jindx[j] + i];
  }

  return e;
}


/* ================================================================
 * Unified incremental MFE recomputation
 * ================================================================ */

PUBLIC float
vrna_mfe_incremental_ex(vrna_fold_compound_t *fc,
                        unsigned int         pos,
                        char                 new_base,
                        unsigned int         flags,
                        char                 *structure,
                        vrna_topk_cache_t    *cache)
{
  unsigned int  i, j, length;
  int           ij, *indx, *f5, *c, *fML, *fM1;
  int           *Fmi, *DMLi;
  float         mfe;
  vrna_md_t     *md;

  int use_lazy   = (flags & VRNA_MFE_INC_LAZY)   ? 1 : 0;
  int use_sparse = (flags & VRNA_MFE_INC_SPARSE)  ? 1 : 0;
  int use_topk   = (flags & VRNA_MFE_INC_TOPK)    ? 1 : 0;

  mfe = (float)(INF / 100.);

  if (!fc || !fc->matrices || !fc->matrices->f5 || !fc->matrices->c) {
    vrna_log_warning("vrna_mfe_incremental_ex: "
                     "fold compound not initialized or missing DP matrices");
    return mfe;
  }

  length = fc->length;

  if (pos < 1 || pos > length) {
    vrna_log_warning("vrna_mfe_incremental_ex: "
                     "mutation position %u out of range [1, %u]",
                     pos, length);
    return mfe;
  }

  if (use_topk && !cache) {
    vrna_log_warning("vrna_mfe_incremental_ex: "
                     "TOPK flag set but no cache provided, falling back");
    use_topk = 0;
  }

  md = &(fc->params->model_details);

  if (md->circ || fc->strands > 1) {
    update_sequence_at_pos(fc, pos, new_base);
    vrna_hc_init(fc);
    return vrna_mfe(fc, structure);
  }

  /* ---- Step 1: Setup (sequence + ptype + hc) ---- */
  if (use_lazy) {
    update_sequence_at_pos_lazy(fc, pos, new_base);
    update_ptype_at_pos(fc, pos);
    update_hc_at_pos(fc, pos);
  } else {
    update_sequence_at_pos(fc, pos, new_base);
    vrna_hc_init(fc);
  }

  /* ---- Step 2: Determine bounds for fML ---- */
  indx = fc->jindx;
  c    = fc->matrices->c;
  fML  = fc->matrices->fML;
  fM1  = fc->matrices->fM1;
  f5   = fc->matrices->f5;

  unsigned int affected_i_max = (pos < length) ? pos + 1 : pos;
  unsigned int fml_i_lo = 1;
  unsigned int fml_i_hi = affected_i_max;

  /* ---- Step 3: Unified DP fill ---- */
  Fmi  = (int *)vrna_alloc(sizeof(int) * (length + 2));
  DMLi = (int *)vrna_alloc(sizeof(int) * (length + 2));

  int change_propagated = 0;

  unsigned int fml_col_start = (pos > 1) ? pos - 1 : 1;

  for (i = length - 1; i >= 1; i--) {
    for (j = 0; j <= length + 1; j++)
      Fmi[j] = DMLi[j] = INF;

    int row_changed = 0;

    for (j = i + 1; j <= length; j++) {
      ij = indx[j] + i;

      if (i <= pos && pos <= j) {
        int old_val = c[ij];
        if (use_topk)
          c[ij] = incremental_decompose_pair_topk(fc, i, j, cache);
        else
          c[ij] = incremental_decompose_pair(fc, i, j);
        if (c[ij] != old_val)
          row_changed = 1;
      }

      int in_spatial_range = (i >= fml_i_lo && i <= fml_i_hi);

      int need_fml = 0;
      if (in_spatial_range) {
        if (use_sparse) {
          if (i == pos || i == pos + 1)
            need_fml = 1;
          else if (i < pos && change_propagated)
            need_fml = 1;
        } else {
          need_fml = 1;
        }

        if (need_fml && j < fml_col_start)
          need_fml = 0;
      }

      if (need_fml) {
        int old_val = fML[ij];
        if (use_topk)
          fML[ij] = topk_E_ml_stems_fast(fc, i, j, Fmi, DMLi, cache);
        else
          fML[ij] = vrna_E_ml_stems_fast(fc, i, j, Fmi, DMLi);
        if (fML[ij] != old_val)
          row_changed = 1;

        if (fM1)
          fM1[ij] = vrna_mfe_multibranch_m1(fc, i, j);
      } else {
        Fmi[j] = fML[ij];
      }
    }

    if (row_changed)
      change_propagated = 1;
  }

  free(Fmi);
  free(DMLi);

  /* ---- Step 4: Recompute f5 ---- */
  (void)vrna_mfe_exterior_f5(fc);

  /* ---- Step 5: Backtrack ---- */
  if (structure) {
    mfe = vrna_backtrack5(fc, length, structure);
    if (mfe >= (float)(INF / 200.)) {
      memset(structure, '.', length);
      structure[length] = '\0';
      mfe = (float)f5[length] / 100.;
    }
  } else {
    mfe = (float)f5[length] / 100.;
  }

  return mfe;
}


PUBLIC float
vrna_mfe_incremental(vrna_fold_compound_t *fc,
                     unsigned int         pos,
                     char                 new_base,
                     char                 *structure)
{
  return vrna_mfe_incremental_ex(fc, pos, new_base, 0, structure, NULL);
}


/* ================================================================
 * Top-K bifurcation heuristic
 * ================================================================ */

static INLINE void
topk_insert(unsigned short *splits,
            int            *energies,
            unsigned char  *count,
            unsigned int   K,
            unsigned short k_val,
            int            energy)
{
  unsigned int c = *count;

  if (c < K) {
    splits[c]   = k_val;
    energies[c] = energy;
    (*count)++;
    unsigned int pos = c;
    while (pos > 0 && energies[pos] > energies[pos - 1]) {
      unsigned short ts = splits[pos]; splits[pos] = splits[pos-1]; splits[pos-1] = ts;
      int te = energies[pos]; energies[pos] = energies[pos-1]; energies[pos-1] = te;
      pos--;
    }
  } else if (energy < energies[0]) {
    splits[0]   = k_val;
    energies[0] = energy;
    unsigned int pos = 0;
    while (pos + 1 < K && energies[pos] < energies[pos + 1]) {
      unsigned short ts = splits[pos]; splits[pos] = splits[pos+1]; splits[pos+1] = ts;
      int te = energies[pos]; energies[pos] = energies[pos+1]; energies[pos+1] = te;
      pos++;
    }
  }
}


static INLINE void
topk_insert_pair(unsigned short *arr_k,
                 unsigned short *arr_l,
                 int            *energies,
                 unsigned char  *count,
                 unsigned int   K,
                 unsigned short k_val,
                 unsigned short l_val,
                 int            energy)
{
  unsigned int c = *count;

  if (c < K) {
    arr_k[c]    = k_val;
    arr_l[c]    = l_val;
    energies[c] = energy;
    (*count)++;
    unsigned int pos = c;
    while (pos > 0 && energies[pos] > energies[pos - 1]) {
      unsigned short tk = arr_k[pos]; arr_k[pos] = arr_k[pos-1]; arr_k[pos-1] = tk;
      unsigned short tl = arr_l[pos]; arr_l[pos] = arr_l[pos-1]; arr_l[pos-1] = tl;
      int te = energies[pos]; energies[pos] = energies[pos-1]; energies[pos-1] = te;
      pos--;
    }
  } else if (energy < energies[0]) {
    arr_k[0]    = k_val;
    arr_l[0]    = l_val;
    energies[0] = energy;
    unsigned int pos = 0;
    while (pos + 1 < K && energies[pos] < energies[pos + 1]) {
      unsigned short tk = arr_k[pos]; arr_k[pos] = arr_k[pos+1]; arr_k[pos+1] = tk;
      unsigned short tl = arr_l[pos]; arr_l[pos] = arr_l[pos+1]; arr_l[pos+1] = tl;
      int te = energies[pos]; energies[pos] = energies[pos+1]; energies[pos+1] = te;
      pos++;
    }
  }
}


PUBLIC vrna_topk_cache_t *
vrna_topk_cache_from_fold(vrna_fold_compound_t *fc,
                          unsigned int         K)
{
  unsigned int      i, j, k, n, ij;
  int               *indx, *fML, e;
  int               *tmp_energies;
  vrna_topk_cache_t *cache;

  if (!fc || !fc->matrices || !fc->matrices->fML) {
    vrna_log_warning("vrna_topk_cache_from_fold: "
                     "fold compound not initialized or missing DP matrices");
    return NULL;
  }

  n    = fc->length;
  indx = fc->jindx;
  fML  = fc->matrices->fML;

  if (K < 1)
    K = 1;
  if (K > 255)
    K = 255;

  cache = (vrna_topk_cache_t *)vrna_alloc(sizeof(vrna_topk_cache_t));
  cache->K    = K;
  cache->n    = n;
  cache->indx = indx;

  size_t matrix_size = (size_t)(n + 1) * (size_t)(n + 2) / 2 + 2;

  cache->fml_splits = (unsigned short *)vrna_alloc(sizeof(unsigned short) * matrix_size * K);
  cache->mb_splits  = (unsigned short *)vrna_alloc(sizeof(unsigned short) * matrix_size * K);
  cache->int_k      = (unsigned short *)vrna_alloc(sizeof(unsigned short) * matrix_size * K);
  cache->int_l      = (unsigned short *)vrna_alloc(sizeof(unsigned short) * matrix_size * K);
  cache->fml_count  = (unsigned char *)calloc(matrix_size, sizeof(unsigned char));
  cache->mb_count   = (unsigned char *)calloc(matrix_size, sizeof(unsigned char));
  cache->int_count  = (unsigned char *)calloc(matrix_size, sizeof(unsigned char));

  tmp_energies = (int *)vrna_alloc(sizeof(int) * K);

  for (i = 1; i < n; i++) {
    for (j = i + 3; j <= n; j++) {
      ij = indx[j] + i;

      cache->fml_count[ij] = 0;
      for (unsigned int t = 0; t < K; t++)
        tmp_energies[t] = INF;

      for (k = i + 1; k <= j - 2; k++) {
        int left  = fML[indx[k] + i];
        int right = fML[indx[j] + k + 1];
        e = ADD_OR_INF(left, right);
        if (e < INF)
          topk_insert(cache->fml_splits + (size_t)ij * K,
                      tmp_energies,
                      &cache->fml_count[ij],
                      K, (unsigned short)k, e);
      }

      cache->mb_count[ij] = 0;
      for (unsigned int t = 0; t < K; t++)
        tmp_energies[t] = INF;

      for (k = i + 2; k <= j - 2; k++) {
        int left  = fML[indx[k] + i + 1];
        int right = fML[indx[j - 1] + k + 1];
        e = ADD_OR_INF(left, right);
        if (e < INF)
          topk_insert(cache->mb_splits + (size_t)ij * K,
                      tmp_energies,
                      &cache->mb_count[ij],
                      K, (unsigned short)k, e);
      }
    }
  }

  /* Interior loop top-K candidates */
  {
    int           *c_arr = fc->matrices->c;
    char          *ptype = fc->ptype;
    short         *S     = fc->sequence_encoding;
    vrna_param_t  *P     = fc->params;
    vrna_md_t     *md    = &(P->model_details);
    vrna_hc_t     *hc    = fc->hc;
    unsigned int  *hc_up = hc->up_int;
    int           noGUclosure = md->noGUclosure;
    int           rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};

    for (i = 1; i < n; i++) {
      for (j = i + 4; j <= n; j++) {
        ij = indx[j] + i;

        unsigned int type_ij = vrna_get_ptype(ij, ptype);
        if (type_ij == 0)
          continue;

        if (noGUclosure && (type_ij == 3 || type_ij == 4))
          continue;

        cache->int_count[ij] = 0;
        for (unsigned int t = 0; t < K; t++)
          tmp_energies[t] = INF;

        unsigned int first_l = i + 2 + 1;
        if (first_l + MAXLOOP + 1 < j)
          first_l = j - 1 - MAXLOOP;

        unsigned int u2 = 1;
        unsigned int l_pos;
        for (l_pos = j - 2; l_pos >= first_l; l_pos--, u2++) {
          if (u2 > hc_up[l_pos + 1])
            break;

          unsigned int last_k = l_pos - 1;
          if (last_k + u2 > i + 1 + MAXLOOP)
            last_k = i + 1 + MAXLOOP - u2;
          if (last_k > i + 1 + hc_up[i + 1])
            last_k = i + 1 + hc_up[i + 1];

          unsigned int u1 = 1;
          unsigned int k_pos;
          for (k_pos = i + 2; k_pos <= last_k; k_pos++, u1++) {
            if (!hc->eval_int(i, j, k_pos, l_pos, hc))
              continue;

            unsigned int kl = indx[l_pos] + k_pos;
            int eee = c_arr[kl];
            if (eee >= INF)
              continue;

            unsigned int type2 = rtype[vrna_get_ptype(kl, ptype)];
            if (noGUclosure && (type2 == 3 || type2 == 4))
              continue;

            eee += vrna_E_internal(u1, u2, type_ij, type2,
                                    S[i + 1], S[j - 1],
                                    S[k_pos - 1], S[l_pos + 1], P);

            if (eee < INF)
              topk_insert_pair(cache->int_k + (size_t)ij * K,
                               cache->int_l + (size_t)ij * K,
                               tmp_energies,
                               &cache->int_count[ij],
                               K,
                               (unsigned short)k_pos,
                               (unsigned short)l_pos,
                               eee);
          }
        }
      }
    }
  }

  free(tmp_energies);
  return cache;
}


PUBLIC void
vrna_topk_cache_free(vrna_topk_cache_t *cache)
{
  if (cache) {
    free(cache->fml_splits);
    free(cache->mb_splits);
    free(cache->int_k);
    free(cache->int_l);
    free(cache->fml_count);
    free(cache->mb_count);
    free(cache->int_count);
    free(cache);
  }
}


static INLINE int
mfe_internal_topk(vrna_fold_compound_t *fc,
                  unsigned int         i,
                  unsigned int         j,
                  vrna_topk_cache_t    *cache)
{
  unsigned int  ij, t, K, kl;
  int           e, eee;
  int           *indx, *c_arr;
  char          *ptype;
  short         *S;
  vrna_param_t  *P;
  vrna_md_t     *md;
  int           rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};

  indx   = fc->jindx;
  c_arr  = fc->matrices->c;
  ptype  = fc->ptype;
  S      = fc->sequence_encoding;
  P      = fc->params;
  md     = &(P->model_details);
  K      = cache->K;
  ij     = indx[j] + i;
  e      = INF;

  unsigned int type_ij = vrna_get_ptype(ij, ptype);

  /* 1. Stack: always evaluate (i+1, j-1) */
  if (j > i + 3) {
    kl = indx[j - 1] + i + 1;
    eee = c_arr[kl];
    if (eee < INF) {
      unsigned int type2 = rtype[vrna_get_ptype(kl, ptype)];
      if (!(md->noGUclosure && (type2 == 3 || type2 == 4))) {
        eee += vrna_E_internal(0, 0, type_ij, type2,
                                S[i + 1], S[j - 1],
                                S[i], S[j], P);
        e = MIN2(e, eee);
      }
    }
  }

  /* 2. Bulges */
  {
    unsigned int max_u, u_pos;
    max_u = (j - i - 3 < MAXLOOP) ? j - i - 3 : MAXLOOP;
    for (u_pos = 1; u_pos <= max_u; u_pos++) {
      unsigned int l_pos = j - 1 - u_pos;
      if (l_pos <= i + 1) break;
      kl = indx[l_pos] + i + 1;
      eee = c_arr[kl];
      if (eee < INF) {
        unsigned int type2 = rtype[vrna_get_ptype(kl, ptype)];
        if (!(md->noGUclosure && (type2 == 3 || type2 == 4))) {
          eee += vrna_E_internal(0, u_pos, type_ij, type2,
                                  S[i + 1], S[j - 1],
                                  S[i], S[l_pos + 1], P);
          e = MIN2(e, eee);
        }
      }
    }
    for (u_pos = 1; u_pos <= max_u; u_pos++) {
      unsigned int k_pos = i + 1 + u_pos;
      if (k_pos >= j - 1) break;
      kl = indx[j - 1] + k_pos;
      eee = c_arr[kl];
      if (eee < INF) {
        unsigned int type2 = rtype[vrna_get_ptype(kl, ptype)];
        if (!(md->noGUclosure && (type2 == 3 || type2 == 4))) {
          eee += vrna_E_internal(u_pos, 0, type_ij, type2,
                                  S[i + 1], S[j - 1],
                                  S[k_pos - 1], S[j], P);
          e = MIN2(e, eee);
        }
      }
    }
  }

  /* 3. Cached general interior loops */
  for (t = 0; t < cache->int_count[ij]; t++) {
    unsigned int k_pos = cache->int_k[(size_t)ij * K + t];
    unsigned int l_pos = cache->int_l[(size_t)ij * K + t];

    kl = indx[l_pos] + k_pos;
    eee = c_arr[kl];
    if (eee >= INF)
      continue;

    unsigned int u1    = k_pos - i - 1;
    unsigned int u2    = j - l_pos - 1;
    unsigned int type2 = rtype[vrna_get_ptype(kl, ptype)];

    if (md->noGUclosure && (type2 == 3 || type2 == 4))
      continue;

    eee += vrna_E_internal(u1, u2, type_ij, type2,
                            S[i + 1], S[j - 1],
                            S[k_pos - 1], S[l_pos + 1], P);
    e = MIN2(e, eee);
  }

  return e;
}


static INLINE int
mb_loop_from_fml_topk(vrna_fold_compound_t *fc,
                      unsigned int         i,
                      unsigned int         j,
                      vrna_topk_cache_t    *cache)
{
  short         *S2;
  unsigned int  tt, n, ij, t, K;
  int           decomp, e, *fML, *indx, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n             = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  S2            = fc->sequence_encoding2;
  indx          = fc->jindx;
  fML           = fc->matrices->fML;
  hc            = fc->hc;

  decomp = INF;

  if (!(hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP))
    return INF;

  if (hc->eval_mb &&
      !hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc))
    return INF;

  ij = indx[j] + i;
  K  = cache->K;

  for (t = 0; t < cache->mb_count[ij]; t++) {
    unsigned int k = cache->mb_splits[(size_t)ij * K + t];
    if (k < i + 2 || k > j - 2)
      continue;
    int fml_left  = fML[indx[k] + i + 1];
    int fml_right = fML[indx[j - 1] + k + 1];
    e = ADD_OR_INF(fml_left, fml_right);
    decomp = MIN2(decomp, e);
  }

  if (decomp != INF) {
    tt = vrna_get_ptype_md(S2[j], S2[i], md);

    if (md->noGUclosure && ((tt == 3) || (tt == 4)))
      return INF;

    switch (dangle_model) {
      case 2:
        decomp += vrna_E_multibranch_stem(tt, S2[j - 1], S2[i + 1], P) +
                  P->MLclosing;
        break;
      default:
        decomp += vrna_E_multibranch_stem(tt, -1, -1, P) +
                  P->MLclosing;
        break;
    }
  }

  return decomp;
}


static INLINE int
incremental_decompose_pair_topk(vrna_fold_compound_t *fc,
                                int                  i,
                                int                  j,
                                vrna_topk_cache_t    *cache)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, dangle_model;

  n             = fc->length;
  dangle_model  = fc->params->model_details.dangles;
  hc_decompose  = fc->hc->mx[n * i + j];
  e             = INF;

  if (hc_decompose) {
    new_c = INF;

    energy  = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
    new_c   = MIN2(new_c, energy);

    energy  = mb_loop_from_fml_topk(fc, i, j, cache);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      energy  = vrna_mfe_multibranch_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    if (cache->int_count && cache->int_k)
      energy = mfe_internal_topk(fc, i, j, cache);
    else
      energy = vrna_mfe_internal(fc, i, j);
    new_c   = MIN2(new_c, energy);

    e = new_c;

    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->c); c++) {
        if (fc->aux_grammar->c[c].cb) {
          energy  = fc->aux_grammar->c[c].cb(fc, i, j, fc->aux_grammar->c[c].data);
          e       = MIN2(e, energy);
        }
      }
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (e != INF))
      e -= fc->pscore[fc->jindx[j] + i];
  }

  return e;
}


static INLINE int
topk_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     int                  *fmi,
                     int                  *dmli,
                     vrna_topk_cache_t    *cache)
{
  short         *S;
  unsigned int  t, K, ij, n, type, dangle_model;
  int           decomp, e, en, *indx, *fm, *c_arr;
  char          *ptype;
  vrna_param_t  *P;
  vrna_hc_t     *hc;

  n             = fc->length;
  indx          = fc->jindx;
  fm            = fc->matrices->fML;
  c_arr         = fc->matrices->c;
  P             = fc->params;
  dangle_model  = P->model_details.dangles;
  S             = fc->sequence_encoding;
  ptype         = fc->ptype;
  hc            = fc->hc;
  K             = cache->K;
  ij            = indx[j] + i;
  e             = INF;

  if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
    en = c_arr[ij];
    if (en != INF) {
      type = vrna_get_ptype(ij, ptype);
      if (dangle_model == 2)
        en += vrna_E_multibranch_stem(type,
                                      (i == 1) ? S[n] : S[i - 1],
                                      S[j + 1], P);
      else
        en += vrna_E_multibranch_stem(type, -1, -1, P);

      e = MIN2(e, en);
    }
  }

  if (i + 1 < j) {
    if (hc->eval_mb &&
        hc->eval_mb(i, j, i + 1, j, VRNA_DECOMP_ML_ML, hc)) {
      en = fm[ij + 1];
      if (en != INF) {
        en += P->MLbase;
        e = MIN2(e, en);
      }
    }
  }

  if (j > i + 1 && fmi[j - 1] != INF) {
    en = fmi[j - 1] + P->MLbase;
    e = MIN2(e, en);
  }

  decomp = INF;
  for (t = 0; t < cache->fml_count[ij]; t++) {
    unsigned int k = cache->fml_splits[(size_t)ij * K + t];
    if (k < i + 1 || k > j - 2)
      continue;
    int left  = fmi[k];
    int right = fm[indx[j] + k + 1];
    en = ADD_OR_INF(left, right);
    decomp = MIN2(decomp, en);
  }

  dmli[j] = decomp;
  e = MIN2(e, decomp);

  fmi[j] = e;

  return e;
}
