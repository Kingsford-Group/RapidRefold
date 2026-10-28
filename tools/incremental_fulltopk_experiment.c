/*
 *  incremental_fulltopk_experiment — Full Top-K heuristic evaluation
 *
 *  Compares exact incremental MFE (ground truth) against the full top-K
 *  heuristic that caches both bifurcation splits AND interior loop pairs.
 *
 *  For each single-nucleotide mutation:
 *    1. Restore DP snapshot from initial fold
 *    2. Run base incremental (exact) as ground truth (gt_mfe)
 *    3. Run full top-K with K=1, K=5, K=10
 *
 *  Output CSV includes MFE, timing, energy difference, and structural
 *  accuracy metrics (sensitivity, PPV, F1) comparing each topk structure
 *  against the ground-truth structure.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/structures/pairtable.h>

#include "mfe_incremental.h"

#define MAX_SEQ 100000
#define INIT_FOLD_REPS 3

/* ---- Timing ---- */
static double
now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}


/* ---- Snapshot ---- */
typedef struct {
  unsigned int  n;
  unsigned int  tri_size;
  unsigned int  lin_size;
  unsigned int  ptype_size;
  unsigned int  hc_size;
  int           *c, *fML, *fM1, *f5;
  short         *seq_enc, *seq_enc2, *enc5, *enc3;
  short         *nuc_enc, *nuc_enc5, *nuc_enc3;
  char          *ptype;
  unsigned char *hc_mx;
  char          *sequence;
  char          *nuc_string;
  char          *structure;
} dp_snapshot_t;


static dp_snapshot_t *
snapshot_alloc(unsigned int n)
{
  dp_snapshot_t *s = (dp_snapshot_t *)calloc(1, sizeof(dp_snapshot_t));
  s->n          = n;
  s->tri_size   = ((n + 1) * (n + 2)) / 2;
  s->lin_size   = n + 2;
  s->ptype_size = (n * (n + 1)) / 2 + 2;
  s->hc_size    = (n + 1) * (n + 1) + 1;

  s->c   = (int *)malloc(sizeof(int) * s->tri_size);
  s->fML = (int *)malloc(sizeof(int) * s->tri_size);
  s->fM1 = (int *)malloc(sizeof(int) * s->tri_size);
  s->f5  = (int *)malloc(sizeof(int) * s->lin_size);

  s->seq_enc  = (short *)malloc(sizeof(short) * s->lin_size);
  s->seq_enc2 = (short *)malloc(sizeof(short) * s->lin_size);
  s->enc5     = (short *)malloc(sizeof(short) * s->lin_size);
  s->enc3     = (short *)malloc(sizeof(short) * s->lin_size);
  s->nuc_enc  = (short *)malloc(sizeof(short) * s->lin_size);
  s->nuc_enc5 = (short *)malloc(sizeof(short) * (n + 1));
  s->nuc_enc3 = (short *)malloc(sizeof(short) * (n + 1));

  s->ptype     = (char *)malloc(sizeof(char) * s->ptype_size);
  s->hc_mx     = (unsigned char *)malloc(sizeof(unsigned char) * s->hc_size);
  s->sequence  = (char *)malloc(n + 1);
  s->nuc_string = (char *)malloc(n + 1);
  s->structure = (char *)malloc(n + 1);

  return s;
}


static void
snapshot_save(dp_snapshot_t *snap, vrna_fold_compound_t *fc, const char *structure)
{
  unsigned int tri = snap->tri_size;
  unsigned int lin = snap->lin_size;

  if (fc->matrices) {
    if (fc->matrices->c)   memcpy(snap->c,   fc->matrices->c,   sizeof(int) * tri);
    if (fc->matrices->fML) memcpy(snap->fML, fc->matrices->fML, sizeof(int) * tri);
    if (fc->matrices->fM1) memcpy(snap->fM1, fc->matrices->fM1, sizeof(int) * tri);
    if (fc->matrices->f5)  memcpy(snap->f5,  fc->matrices->f5,  sizeof(int) * lin);
  }

  if (fc->sequence_encoding)  memcpy(snap->seq_enc,  fc->sequence_encoding,  sizeof(short) * lin);
  if (fc->sequence_encoding2) memcpy(snap->seq_enc2, fc->sequence_encoding2, sizeof(short) * lin);
  if (fc->encoding5)          memcpy(snap->enc5,     fc->encoding5,          sizeof(short) * lin);
  if (fc->encoding3)          memcpy(snap->enc3,     fc->encoding3,          sizeof(short) * lin);

  if (fc->ptype)              memcpy(snap->ptype, fc->ptype, sizeof(char) * snap->ptype_size);
  if (fc->hc && fc->hc->mx)  memcpy(snap->hc_mx, fc->hc->mx, sizeof(unsigned char) * snap->hc_size);
  if (fc->sequence)           memcpy(snap->sequence, fc->sequence, snap->n + 1);

  if (fc->nucleotides) {
    if (fc->nucleotides[0].encoding)
      memcpy(snap->nuc_enc, fc->nucleotides[0].encoding, sizeof(short) * lin);
    if (fc->nucleotides[0].encoding5)
      memcpy(snap->nuc_enc5, fc->nucleotides[0].encoding5, sizeof(short) * (snap->n + 1));
    if (fc->nucleotides[0].encoding3)
      memcpy(snap->nuc_enc3, fc->nucleotides[0].encoding3, sizeof(short) * (snap->n + 1));
    if (fc->nucleotides[0].string)
      memcpy(snap->nuc_string, fc->nucleotides[0].string, snap->n + 1);
  }

  if (structure) memcpy(snap->structure, structure, snap->n + 1);
}


static void
snapshot_restore(dp_snapshot_t *snap, vrna_fold_compound_t *fc, char *structure)
{
  unsigned int tri = snap->tri_size;
  unsigned int lin = snap->lin_size;

  if (fc->matrices) {
    if (fc->matrices->c)   memcpy(fc->matrices->c,   snap->c,   sizeof(int) * tri);
    if (fc->matrices->fML) memcpy(fc->matrices->fML, snap->fML, sizeof(int) * tri);
    if (fc->matrices->fM1) memcpy(fc->matrices->fM1, snap->fM1, sizeof(int) * tri);
    if (fc->matrices->f5)  memcpy(fc->matrices->f5,  snap->f5,  sizeof(int) * lin);
  }

  if (fc->sequence_encoding)  memcpy(fc->sequence_encoding,  snap->seq_enc,  sizeof(short) * lin);
  if (fc->sequence_encoding2) memcpy(fc->sequence_encoding2, snap->seq_enc2, sizeof(short) * lin);
  if (fc->encoding5)          memcpy(fc->encoding5,          snap->enc5,     sizeof(short) * lin);
  if (fc->encoding3)          memcpy(fc->encoding3,          snap->enc3,     sizeof(short) * lin);

  if (fc->ptype)              memcpy(fc->ptype, snap->ptype, sizeof(char) * snap->ptype_size);
  if (fc->hc && fc->hc->mx)  memcpy(fc->hc->mx, snap->hc_mx, sizeof(unsigned char) * snap->hc_size);
  if (fc->sequence)           memcpy(fc->sequence, snap->sequence, snap->n + 1);

  if (fc->nucleotides) {
    if (fc->nucleotides[0].encoding)
      memcpy(fc->nucleotides[0].encoding, snap->nuc_enc, sizeof(short) * lin);
    if (fc->nucleotides[0].encoding5)
      memcpy(fc->nucleotides[0].encoding5, snap->nuc_enc5, sizeof(short) * (snap->n + 1));
    if (fc->nucleotides[0].encoding3)
      memcpy(fc->nucleotides[0].encoding3, snap->nuc_enc3, sizeof(short) * (snap->n + 1));
    if (fc->nucleotides[0].string)
      memcpy(fc->nucleotides[0].string, snap->nuc_string, snap->n + 1);
  }

  if (structure) memcpy(structure, snap->structure, snap->n + 1);
}


static void
snapshot_free(dp_snapshot_t *snap)
{
  if (!snap) return;
  free(snap->c); free(snap->fML); free(snap->fM1); free(snap->f5);
  free(snap->seq_enc); free(snap->seq_enc2);
  free(snap->enc5); free(snap->enc3);
  free(snap->nuc_enc); free(snap->nuc_enc5); free(snap->nuc_enc3);
  free(snap->ptype); free(snap->hc_mx);
  free(snap->sequence); free(snap->nuc_string);
  free(snap->structure);
  free(snap);
}


/* ---- Structural comparison metrics ---- */
typedef struct {
  int     tp, fp, fn;
  double  sensitivity;  /* TP / (TP + FN) */
  double  ppv;          /* TP / (TP + FP) */
  double  f1;           /* 2 * sens * ppv / (sens + ppv) */
} struct_metrics_t;

static struct_metrics_t
compare_structures(const char *ref, const char *pred, unsigned int n)
{
  struct_metrics_t m = {0, 0, 0, 1.0, 1.0, 1.0};

  short *pt_ref  = vrna_ptable(ref);
  short *pt_pred = vrna_ptable(pred);

  for (unsigned int i = 1; i <= n; i++) {
    int ref_paired  = (pt_ref[i] > (short)i);
    int pred_paired = (pt_pred[i] > (short)i);

    if (ref_paired && pred_paired && pt_ref[i] == pt_pred[i])
      m.tp++;
    else {
      if (ref_paired)  m.fn++;
      if (pred_paired) m.fp++;
    }
  }

  free(pt_ref);
  free(pt_pred);

  if (m.tp + m.fn > 0)
    m.sensitivity = (double)m.tp / (m.tp + m.fn);
  else
    m.sensitivity = 1.0;

  if (m.tp + m.fp > 0)
    m.ppv = (double)m.tp / (m.tp + m.fp);
  else
    m.ppv = 1.0;

  if (m.sensitivity + m.ppv > 0.0)
    m.f1 = 2.0 * m.sensitivity * m.ppv / (m.sensitivity + m.ppv);
  else
    m.f1 = 0.0;

  return m;
}


/* ---- FASTA reader ---- */
typedef struct { char *id; char *seq; } fasta_entry_t;

static int
read_fasta(const char *path, fasta_entry_t *entries, int max_entries)
{
  FILE *fp = fopen(path, "r");
  if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 0; }

  char line[MAX_SEQ];
  int n = 0;
  char *cur_id = NULL;
  char *cur_seq = NULL;
  size_t seq_len = 0;

  while (fgets(line, sizeof(line), fp)) {
    size_t len = strlen(line);
    while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
      line[--len] = '\0';

    if (line[0] == '>') {
      if (cur_id && cur_seq && n < max_entries) {
        entries[n].id = cur_id;
        entries[n].seq = cur_seq;
        n++;
      } else {
        free(cur_id); free(cur_seq);
      }
      cur_id = strdup(line + 1);
      cur_seq = (char *)calloc(MAX_SEQ, 1);
      seq_len = 0;
    } else if (cur_seq) {
      memcpy(cur_seq + seq_len, line, len);
      seq_len += len;
      cur_seq[seq_len] = '\0';
    }
  }

  if (cur_id && cur_seq && n < max_entries) {
    entries[n].id = cur_id;
    entries[n].seq = cur_seq;
    n++;
  } else {
    free(cur_id); free(cur_seq);
  }

  fclose(fp);
  return n;
}


static void
print_help(const char *progname)
{
  fprintf(stderr,
    "Usage: %s input.fa [options]\n"
    "\n"
    "Saturation mutagenesis experiment comparing exact incremental MFE\n"
    "(ground truth) against the full top-K heuristic (bifurcation +\n"
    "interior loop caching) at K=1,5,10.\n"
    "\n"
    "For each single-nucleotide mutation, reports MFE energy, timing, energy\n"
    "difference, and structural accuracy (sensitivity, PPV, F1) of the top-K\n"
    "predicted structure vs the ground-truth structure.\n"
    "\n"
    "Options:\n"
    "  -o FILE   Write CSV output to FILE (default: stdout)\n"
    "  -s PCT    Sample PCT%% of positions (1-100, default: 100 = all)\n"
    "            Positions are evenly spaced (stratified) for balanced coverage.\n"
    "  -h        Show this help message and exit\n"
    "\n"
    "Input:  FASTA file with one or more RNA sequences.\n"
    "Output: CSV with columns:\n"
    "  seq_id, length, pos, orig, new, method, mfe, time_ms,\n"
    "  mfe_diff, sensitivity, ppv, f1\n"
    "\n"
    "Methods in output:\n"
    "  initial_fold      Initial MFE fold (median of %d repeats)\n"
    "  topk_cache_KN     Time to build full top-K cache for K=N\n"
    "  gt_mfe            Ground-truth incremental MFE (exact)\n"
    "  topk_full_KN      Full top-K heuristic with K=N\n"
    "\n"
    "Examples:\n"
    "  %s sequences.fa\n"
    "  %s sequences.fa -o results.csv\n"
    "  %s sequences.fa -s 10           # 10%% stratified sample\n",
    progname, INIT_FOLD_REPS, progname, progname, progname);
}


int
main(int argc, char *argv[])
{
  const char *infile = NULL;
  const char *outfile = NULL;
  int sample_pct = 100;
  int opt;

  while ((opt = getopt(argc, argv, "ho:s:")) != -1) {
    switch (opt) {
      case 'h':
        print_help(argv[0]);
        return 0;
      case 'o':
        outfile = optarg;
        break;
      case 's':
        sample_pct = atoi(optarg);
        if (sample_pct < 1 || sample_pct > 100) {
          fprintf(stderr, "Error: -s must be between 1 and 100.\n");
          return 1;
        }
        break;
      default:
        print_help(argv[0]);
        return 1;
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: no input file specified.\n\n");
    print_help(argv[0]);
    return 1;
  }

  infile = argv[optind];

  fasta_entry_t entries[1000];
  int n_entries = read_fasta(infile, entries, 1000);
  fprintf(stderr, "Read %d sequences from %s\n", n_entries, infile);

  FILE *out = outfile ? fopen(outfile, "w") : stdout;
  if (!out) { fprintf(stderr, "Cannot open %s\n", outfile); return 1; }

  /* CSV header */
  fprintf(out, "seq_id,length,pos,orig,new,method,mfe,time_ms,mfe_diff,sensitivity,ppv,f1\n");

  const char bases[] = "ACGU";
  const int K_values[] = {1, 5, 10};
  const int n_K = 3;

  for (int si = 0; si < n_entries; si++) {
    char *seq_id = entries[si].id;
    char *seq    = entries[si].seq;
    unsigned int seqlen = strlen(seq);

    /* Build list of positions to test */
    unsigned int n_positions;
    unsigned int *positions;

    if (sample_pct < 100) {
      n_positions = (unsigned int)ceil(seqlen * sample_pct / 100.0);
      if (n_positions < 1) n_positions = 1;
      if (n_positions > seqlen) n_positions = seqlen;
      positions = (unsigned int *)malloc(sizeof(unsigned int) * n_positions);
      for (unsigned int i = 0; i < n_positions; i++) {
        /* Stratified: pick center of each equal-sized bucket */
        positions[i] = (unsigned int)floor((i + 0.5) * seqlen / n_positions) + 1;
        if (positions[i] > seqlen) positions[i] = seqlen;
      }
    } else {
      n_positions = seqlen;
      positions = (unsigned int *)malloc(sizeof(unsigned int) * n_positions);
      for (unsigned int i = 0; i < n_positions; i++)
        positions[i] = i + 1;
    }

    fprintf(stderr, "Processing %s (n=%u, %u positions, %u mutations%s)...\n",
            seq_id, seqlen, n_positions, n_positions * 3,
            sample_pct < 100 ? ", stratified" : "");

    /* ---- Initial fold with INIT_FOLD_REPS repeats, median timing ---- */
    double init_times[INIT_FOLD_REPS];
    char *init_struct = (char *)vrna_alloc(seqlen + 1);
    float init_mfe = 0;

    for (int rep = 0; rep < INIT_FOLD_REPS; rep++) {
      vrna_fold_compound_t *fc_rep = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
      char *s_rep = (char *)vrna_alloc(seqlen + 1);

      double t0 = now_sec();
      float mfe_rep = vrna_mfe(fc_rep, s_rep);
      init_times[rep] = (now_sec() - t0) * 1e3;

      if (rep == 0) {
        init_mfe = mfe_rep;
        memcpy(init_struct, s_rep, seqlen + 1);
      }

      vrna_fold_compound_free(fc_rep);
      free(s_rep);
    }

    /* Bubble sort to find median */
    for (int i = 0; i < INIT_FOLD_REPS - 1; i++)
      for (int j = i + 1; j < INIT_FOLD_REPS; j++) {
        if (init_times[j] < init_times[i]) {
          double tmp = init_times[i];
          init_times[i] = init_times[j];
          init_times[j] = tmp;
        }
      }
    double t_init = init_times[INIT_FOLD_REPS / 2];

    fprintf(stderr, "  Initial fold: %.2f kcal/mol in %.1f ms (median of %d)\n",
            init_mfe, t_init, INIT_FOLD_REPS);

    /* Write initial fold row */
    fprintf(out, "%s,%u,0,-,-,initial_fold,%.2f,%.1f,0.00,1.000,1.000,1.000\n",
            seq_id, seqlen, init_mfe, t_init);

    /* Create the fold compound we'll reuse */
    vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
    init_mfe = vrna_mfe(fc, init_struct);

    /* Save snapshot */
    dp_snapshot_t *snap = snapshot_alloc(seqlen);
    snapshot_save(snap, fc, init_struct);

    /* Build full top-K caches (bifurcation + interior loops) for each K */
    vrna_topk_cache_t *caches[3];
    for (int ki = 0; ki < n_K; ki++) {
      snapshot_restore(snap, fc, init_struct);
      double t0 = now_sec();
      caches[ki] = vrna_topk_cache_from_fold(fc, K_values[ki]);
      double t_cache = (now_sec() - t0) * 1e3;
      fprintf(stderr, "  Full cache K=%d built in %.1f ms\n",
              K_values[ki], t_cache);
      fprintf(out, "%s,%u,0,-,-,topk_cache_K%d,0,%.1f,0.00,,,\n",
              seq_id, seqlen, K_values[ki], t_cache);
    }

    char *work_struct = (char *)vrna_alloc(seqlen + 1);
    char *gt_struct   = (char *)vrna_alloc(seqlen + 1);
    char *topk_struct = (char *)vrna_alloc(seqlen + 1);

    /* Process mutations at selected positions */
    for (unsigned int pi = 0; pi < n_positions; pi++) {
      unsigned int pos = positions[pi];
      char orig = seq[pos - 1];
      char orig_upper = toupper(orig);

      for (int b = 0; b < 4; b++) {
        if (bases[b] == orig_upper) continue;

        /* ---- gt_mfe: base incremental (exact) ---- */
        snapshot_restore(snap, fc, work_struct);
        double t0 = now_sec();
        float gt_mfe = vrna_mfe_incremental(fc, pos, bases[b], gt_struct);
        double t_gt = (now_sec() - t0) * 1e3;

        fprintf(out, "%s,%u,%u,%c,%c,gt_mfe,%.2f,%.3f,0.00,1.000,1.000,1.000\n",
                seq_id, seqlen, pos, orig, bases[b], gt_mfe, t_gt);

        /* ---- full topk for each K ---- */
        for (int ki = 0; ki < n_K; ki++) {
          snapshot_restore(snap, fc, work_struct);
          t0 = now_sec();
          float e_topk = vrna_mfe_incremental_ex(fc, pos, bases[b],
                                                  VRNA_MFE_INC_TOPK,
                                                  topk_struct,
                                                  caches[ki]);
          double t_topk = (now_sec() - t0) * 1e3;

          struct_metrics_t sm = compare_structures(gt_struct, topk_struct, seqlen);

          fprintf(out, "%s,%u,%u,%c,%c,topk_full_K%d,%.2f,%.3f,%.2f,%.3f,%.3f,%.3f\n",
                  seq_id, seqlen, pos, orig, bases[b],
                  K_values[ki], e_topk, t_topk,
                  fabs(e_topk - gt_mfe),
                  sm.sensitivity, sm.ppv, sm.f1);
        }
      }

      if (n_positions >= 10) {
        fprintf(stderr, "\r  pos %u/%u (%u/%u)", pos, seqlen, pi + 1, n_positions);
        fflush(stderr);
      }
    }
    if (n_positions >= 10) fprintf(stderr, "\n");

    /* Cleanup */
    for (int ki = 0; ki < n_K; ki++)
      vrna_topk_cache_free(caches[ki]);
    snapshot_free(snap);
    vrna_fold_compound_free(fc);
    free(positions);
    free(init_struct);
    free(work_struct);
    free(gt_struct);
    free(topk_struct);
  }

  if (outfile) fclose(out);

  for (int i = 0; i < n_entries; i++) {
    free(entries[i].id);
    free(entries[i].seq);
  }

  fprintf(stderr, "Done.\n");
  return 0;
}
