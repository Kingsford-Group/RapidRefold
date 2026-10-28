/*
 *  Tests for incremental MFE update methods.
 *
 *  Exact methods are verified against full recomputation from scratch.
 *  Heuristic methods (topk) are compared against the exact incremental
 *  method and report accuracy metrics.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/mfe/global.h>

#include "mfe_incremental.h"

#define FLOAT_EPS 0.011

static int total_tests = 0;
static int passed_tests = 0;

static void
test_pass(const char *name)
{
  total_tests++;
  passed_tests++;
  printf("  PASS: %s\n", name);
}

static void
test_fail(const char *name, const char *detail)
{
  total_tests++;
  printf("  FAIL: %s -- %s\n", name, detail);
}


/* ---- Basic exact incremental tests ---- */

static void
test_basic_incremental(void)
{
  printf("\n--- Basic incremental tests ---\n");

  struct {
    const char *seq;
    unsigned int pos;
    char new_base;
    const char *label;
  } cases[] = {
    {"GGGAAACCC",     5, 'U', "Middle mutation"},
    {"GGGAAACCC",     1, 'A', "5' end mutation"},
    {"GGGAAACCC",     9, 'A', "3' end mutation"},
    {"GGGGAAAACCCC", 3, 'A', "Break base pair"},
    {"GGGAAACCC",     5, 'C', "Hairpin loop mutation"},
    {"GCGCGCAAAGCGCGC", 8, 'U', "GC-rich stem mutation"},
  };
  int n_cases = 6;

  for (int t = 0; t < n_cases; t++) {
    unsigned int n = strlen(cases[t].seq);

    /* Full recomp */
    char *seq2 = strdup(cases[t].seq);
    seq2[cases[t].pos - 1] = cases[t].new_base;
    vrna_fold_compound_t *fc2 = vrna_fold_compound(seq2, NULL, VRNA_OPTION_MFE);
    char *s2 = (char *)vrna_alloc(n + 1);
    float e2 = vrna_mfe(fc2, s2);

    /* Incremental */
    vrna_fold_compound_t *fc = vrna_fold_compound(cases[t].seq, NULL, VRNA_OPTION_MFE);
    char *s = (char *)vrna_alloc(n + 1);
    vrna_mfe(fc, s);
    float e = vrna_mfe_incremental(fc, cases[t].pos, cases[t].new_base, s);

    if (fabs(e - e2) < FLOAT_EPS && strcmp(s, s2) == 0)
      test_pass(cases[t].label);
    else {
      char buf[256];
      snprintf(buf, sizeof(buf), "full=%.2f inc=%.2f", e2, e);
      test_fail(cases[t].label, buf);
    }

    vrna_fold_compound_free(fc);
    vrna_fold_compound_free(fc2);
    free(s); free(s2); free(seq2);
  }
}


/* ---- Exhaustive exact method test (100bp) ---- */

static void
test_exact_exhaustive(void)
{
  printf("\n--- Exhaustive 100bp: incremental ---\n");

  srand(42);
  const char bases[] = "ACGU";
  int n = 100;
  char *seq = malloc(n + 1);
  for (int i = 0; i < n; i++) seq[i] = bases[rand() % 4];
  seq[n] = '\0';

  int pass = 0, total = 0;

  for (int pos = 1; pos <= n; pos++) {
    for (int b = 0; b < 4; b++) {
      if (bases[b] == seq[pos - 1]) continue;
      total++;

      char *seq2 = strdup(seq);
      seq2[pos - 1] = bases[b];
      vrna_fold_compound_t *fc2 = vrna_fold_compound(seq2, NULL, VRNA_OPTION_MFE);
      char *s2 = (char *)vrna_alloc(n + 1);
      float e2 = vrna_mfe(fc2, s2);

      vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
      char *s = (char *)vrna_alloc(n + 1);
      vrna_mfe(fc, s);
      float e = vrna_mfe_incremental(fc, pos, bases[b], s);

      if (fabs(e - e2) < FLOAT_EPS && strcmp(s, s2) == 0)
        pass++;

      vrna_fold_compound_free(fc);
      vrna_fold_compound_free(fc2);
      free(s); free(s2); free(seq2);
    }
  }

  char buf[128];
  snprintf(buf, sizeof(buf), "incremental: %d/%d", pass, total);
  if (pass == total)
    test_pass(buf);
  else
    test_fail(buf, "not all mutations matched full recomp");

  free(seq);
}


/* ---- Exhaustive heuristic test: topk (100bp) ---- */

static void
test_topk_exhaustive(void)
{
  printf("\n--- Exhaustive 100bp heuristic: topk5 ---\n");

  srand(42);
  const char bases[] = "ACGU";
  int n = 100;
  char *seq = malloc(n + 1);
  for (int i = 0; i < n; i++) seq[i] = bases[rand() % 4];
  seq[n] = '\0';

  /* Build top-K cache */
  vrna_fold_compound_t *fcc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
  char *sc = (char *)vrna_alloc(n + 1);
  vrna_mfe(fcc, sc);
  vrna_topk_cache_t *cache = vrna_topk_cache_from_fold(fcc, 5);
  vrna_fold_compound_free(fcc);
  free(sc);

  int exact_match = 0, total = 0;

  for (int pos = 1; pos <= n; pos++) {
    for (int b = 0; b < 4; b++) {
      if (bases[b] == seq[pos - 1]) continue;
      total++;

      /* Reference: exact incremental */
      vrna_fold_compound_t *fc_ref = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
      char *s_ref = (char *)vrna_alloc(n + 1);
      vrna_mfe(fc_ref, s_ref);
      float e_ref = vrna_mfe_incremental(fc_ref, pos, bases[b], s_ref);

      /* TopK */
      vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
      char *s = (char *)vrna_alloc(n + 1);
      vrna_mfe(fc, s);
      float e_test = vrna_mfe_incremental_ex(fc, pos, bases[b],
                                              VRNA_MFE_INC_TOPK, s, cache);

      if (fabs(e_test - e_ref) < FLOAT_EPS && strcmp(s, s_ref) == 0)
        exact_match++;

      vrna_fold_compound_free(fc);
      vrna_fold_compound_free(fc_ref);
      free(s); free(s_ref);
    }
  }

  vrna_topk_cache_free(cache);

  char buf[128];
  snprintf(buf, sizeof(buf), "topk5: %d/%d exact matches (%.1f%%)",
           exact_match, total, 100.0 * exact_match / total);
  test_pass(buf);

  free(seq);
}


/* ---- Top-K cache test ---- */

static void
test_topk_cache(void)
{
  printf("\n--- Top-K cache ---\n");

  const char *seq = "GGGAAACCCUUU";
  unsigned int n = strlen(seq);
  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
  char *s = (char *)vrna_alloc(n + 1);
  vrna_mfe(fc, s);

  vrna_topk_cache_t *cache = vrna_topk_cache_from_fold(fc, 5);
  if (cache && cache->K == 5 && cache->n == n) {
    char buf[64];
    snprintf(buf, sizeof(buf), "Cache: K=%u, n=%u", cache->K, cache->n);
    test_pass(buf);
  } else {
    test_fail("Top-K cache creation", "NULL or wrong dimensions");
  }

  vrna_topk_cache_free(cache);
  vrna_fold_compound_free(fc);
  free(s);
}


/* ---- Composable flag combination tests ---- */

static void
test_composable_flags(void)
{
  printf("\n--- Composable flag combinations (100bp exhaustive) ---\n");

  srand(42);
  const char bases_c[] = "ACGU";
  int n_c = 100;
  char *seq_c = malloc(n_c + 1);
  for (int i = 0; i < n_c; i++) seq_c[i] = bases_c[rand() % 4];
  seq_c[n_c] = '\0';

  /* Build top-K cache for topk combos */
  vrna_fold_compound_t *fc_cache = vrna_fold_compound(seq_c, NULL, VRNA_OPTION_MFE);
  char *s_cache = (char *)vrna_alloc(n_c + 1);
  vrna_mfe(fc_cache, s_cache);
  vrna_topk_cache_t *cache_c = vrna_topk_cache_from_fold(fc_cache, 5);
  vrna_fold_compound_free(fc_cache);
  free(s_cache);

  typedef struct {
    const char   *name;
    unsigned int  flags;
    int           exact;
    int           need_cache;
  } combo_t;

  combo_t combos[] = {
    /* Exact combinations */
    {"lazy_only",           VRNA_MFE_INC_LAZY,                                    1, 0},
    {"sparse_only",         VRNA_MFE_INC_SPARSE,                                  1, 0},
    {"lazy+sparse",         VRNA_MFE_INC_LAZY | VRNA_MFE_INC_SPARSE,              1, 0},
    /* Heuristic combinations */
    {"topk",                VRNA_MFE_INC_TOPK,                                    0, 1},
    {"lazy+topk",           VRNA_MFE_INC_LAZY | VRNA_MFE_INC_TOPK,               0, 1},
    {"sparse+topk",         VRNA_MFE_INC_SPARSE | VRNA_MFE_INC_TOPK,             0, 1},
    {"lazy+sparse+topk",    VRNA_MFE_INC_LAZY | VRNA_MFE_INC_SPARSE | VRNA_MFE_INC_TOPK, 0, 1},
  };
  int n_combos = sizeof(combos) / sizeof(combos[0]);

  for (int ci = 0; ci < n_combos; ci++) {
    combo_t *co = &combos[ci];
    int pass = 0, total = 0;

    for (int pos = 1; pos <= n_c; pos++) {
      for (int b = 0; b < 4; b++) {
        if (bases_c[b] == seq_c[pos - 1]) continue;
        total++;

        float e_ref;
        char *s_ref;
        vrna_fold_compound_t *fc_ref;

        if (co->exact) {
          char *seq2 = strdup(seq_c);
          seq2[pos - 1] = bases_c[b];
          fc_ref = vrna_fold_compound(seq2, NULL, VRNA_OPTION_MFE);
          s_ref = (char *)vrna_alloc(n_c + 1);
          e_ref = vrna_mfe(fc_ref, s_ref);
          free(seq2);
        } else {
          fc_ref = vrna_fold_compound(seq_c, NULL, VRNA_OPTION_MFE);
          s_ref = (char *)vrna_alloc(n_c + 1);
          vrna_mfe(fc_ref, s_ref);
          e_ref = vrna_mfe_incremental(fc_ref, pos, bases_c[b], s_ref);
        }

        vrna_fold_compound_t *fc = vrna_fold_compound(seq_c, NULL, VRNA_OPTION_MFE);
        char *s = (char *)vrna_alloc(n_c + 1);
        vrna_mfe(fc, s);

        float e_test = vrna_mfe_incremental_ex(fc, pos, bases_c[b],
                                                co->flags, s,
                                                co->need_cache ? cache_c : NULL);

        if (fabs(e_test - e_ref) < FLOAT_EPS && strcmp(s, s_ref) == 0)
          pass++;

        vrna_fold_compound_free(fc);
        vrna_fold_compound_free(fc_ref);
        free(s); free(s_ref);
      }
    }

    char buf[256];
    if (co->exact) {
      snprintf(buf, sizeof(buf), "%s: %d/%d", co->name, pass, total);
      if (pass == total)
        test_pass(buf);
      else
        test_fail(buf, "not all mutations matched full recomp");
    } else {
      snprintf(buf, sizeof(buf), "%s: %d/%d exact (%.1f%%)",
               co->name, pass, total, 100.0 * pass / total);
      test_pass(buf);
    }
  }

  vrna_topk_cache_free(cache_c);
  free(seq_c);
}


int
main(void)
{
  printf("=== Incremental MFE Update Tests ===\n");

  test_basic_incremental();
  test_exact_exhaustive();
  test_topk_exhaustive();
  test_topk_cache();
  test_composable_flags();

  printf("\n=== Summary ===\n");
  printf("Total: %d, Passed: %d, Failed: %d\n",
         total_tests, passed_tests, total_tests - passed_tests);

  return (passed_tests == total_tests) ? 0 : 1;
}
