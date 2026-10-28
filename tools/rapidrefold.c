/*
 *  RapidRefold — Fast incremental RNA refolding under point mutations.
 *
 *  Two subcommands:
 *    cache  — fold a sequence and serialize DP state to a binary file
 *    fold   — load a cache and refold under specified mutations
 *
 *  Build with: make rapidrefold
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <pthread.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/structures/pairtable.h>

#include "mfe_incremental.h"

#define MAX_SEQ     100000
#define MAGIC       "RRFC"
#define VERSION     1
#define MODE_INC    0
#define MODE_TOPK   1

/* ================================================================
 * Timing
 * ================================================================ */

static double
now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}


/* ================================================================
 * DP snapshot — save/restore fold compound state
 * ================================================================ */

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

  return s;
}


static void
snapshot_save(dp_snapshot_t *snap, vrna_fold_compound_t *fc)
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
}


static void
snapshot_restore(dp_snapshot_t *snap, vrna_fold_compound_t *fc)
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
  free(snap);
}


/* ================================================================
 * Binary cache file I/O
 * ================================================================ */

static int
write_cache_file(const char *path,
                 int mode,
                 const char *sequence,
                 const char *structure,
                 float init_mfe,
                 unsigned int K,
                 dp_snapshot_t *snap,
                 vrna_topk_cache_t *topk_cache)
{
  FILE *fp = fopen(path, "wb");
  if (!fp) { fprintf(stderr, "Error: cannot open %s for writing\n", path); return -1; }

  unsigned int n = snap->n;
  uint32_t ver = VERSION;
  uint8_t  m   = (uint8_t)mode;
  uint32_t len = n;
  uint32_t kval = K;
  int32_t  mfe_int = (int32_t)(init_mfe * 100.0f);

  fwrite(MAGIC, 1, 4, fp);
  fwrite(&ver, sizeof(uint32_t), 1, fp);
  fwrite(&m, sizeof(uint8_t), 1, fp);
  fwrite(&len, sizeof(uint32_t), 1, fp);
  fwrite(&kval, sizeof(uint32_t), 1, fp);
  fwrite(&mfe_int, sizeof(int32_t), 1, fp);

  fwrite(sequence, 1, n, fp);
  fwrite(structure, 1, n, fp);

  /* DP snapshot */
  fwrite(snap->c,        sizeof(int),           snap->tri_size,   fp);
  fwrite(snap->fML,      sizeof(int),           snap->tri_size,   fp);
  fwrite(snap->fM1,      sizeof(int),           snap->tri_size,   fp);
  fwrite(snap->f5,       sizeof(int),           snap->lin_size,   fp);
  fwrite(snap->seq_enc,  sizeof(short),         snap->lin_size,   fp);
  fwrite(snap->seq_enc2, sizeof(short),         snap->lin_size,   fp);
  fwrite(snap->enc5,     sizeof(short),         snap->lin_size,   fp);
  fwrite(snap->enc3,     sizeof(short),         snap->lin_size,   fp);
  fwrite(snap->nuc_enc,  sizeof(short),         snap->lin_size,   fp);
  fwrite(snap->nuc_enc5, sizeof(short),         n + 1,            fp);
  fwrite(snap->nuc_enc3, sizeof(short),         n + 1,            fp);
  fwrite(snap->ptype,    sizeof(char),          snap->ptype_size, fp);
  fwrite(snap->hc_mx,    sizeof(unsigned char), snap->hc_size,    fp);
  fwrite(snap->sequence, 1,                     n + 1,            fp);
  fwrite(snap->nuc_string, 1,                   n + 1,            fp);

  /* Top-K cache arrays (only in topk mode) */
  if (mode == MODE_TOPK && topk_cache) {
    size_t matrix_size = (size_t)(n + 1) * (size_t)(n + 2) / 2 + 2;
    fwrite(topk_cache->fml_splits, sizeof(unsigned short), matrix_size * K, fp);
    fwrite(topk_cache->mb_splits,  sizeof(unsigned short), matrix_size * K, fp);
    fwrite(topk_cache->int_k,      sizeof(unsigned short), matrix_size * K, fp);
    fwrite(topk_cache->int_l,      sizeof(unsigned short), matrix_size * K, fp);
    fwrite(topk_cache->fml_count,  sizeof(unsigned char),  matrix_size,     fp);
    fwrite(topk_cache->mb_count,   sizeof(unsigned char),  matrix_size,     fp);
    fwrite(topk_cache->int_count,  sizeof(unsigned char),  matrix_size,     fp);
  }

  fclose(fp);
  return 0;
}


typedef struct {
  int               mode;
  unsigned int      n;
  unsigned int      K;
  float             init_mfe;
  char              *sequence;
  char              *structure;
  dp_snapshot_t     *snap;
  vrna_topk_cache_t *topk_cache;
} cache_data_t;


static cache_data_t *
read_cache_file(const char *path)
{
  FILE *fp = fopen(path, "rb");
  if (!fp) { fprintf(stderr, "Error: cannot open %s for reading\n", path); return NULL; }

  char magic[4];
  uint32_t ver;
  uint8_t  m;
  uint32_t len, kval;
  int32_t  mfe_int;

  if (fread(magic, 1, 4, fp) != 4 || memcmp(magic, MAGIC, 4) != 0) {
    fprintf(stderr, "Error: %s is not a valid RapidRefold cache file\n", path);
    fclose(fp); return NULL;
  }

  fread(&ver, sizeof(uint32_t), 1, fp);
  if (ver != VERSION) {
    fprintf(stderr, "Error: cache version %u not supported (expected %d)\n", ver, VERSION);
    fclose(fp); return NULL;
  }

  fread(&m, sizeof(uint8_t), 1, fp);
  fread(&len, sizeof(uint32_t), 1, fp);
  fread(&kval, sizeof(uint32_t), 1, fp);
  fread(&mfe_int, sizeof(int32_t), 1, fp);

  unsigned int n = len;
  unsigned int K = kval;

  char *sequence  = (char *)malloc(n + 1);
  char *structure = (char *)malloc(n + 1);
  fread(sequence, 1, n, fp);  sequence[n] = '\0';
  fread(structure, 1, n, fp); structure[n] = '\0';

  dp_snapshot_t *snap = snapshot_alloc(n);

  fread(snap->c,        sizeof(int),           snap->tri_size,   fp);
  fread(snap->fML,      sizeof(int),           snap->tri_size,   fp);
  fread(snap->fM1,      sizeof(int),           snap->tri_size,   fp);
  fread(snap->f5,       sizeof(int),           snap->lin_size,   fp);
  fread(snap->seq_enc,  sizeof(short),         snap->lin_size,   fp);
  fread(snap->seq_enc2, sizeof(short),         snap->lin_size,   fp);
  fread(snap->enc5,     sizeof(short),         snap->lin_size,   fp);
  fread(snap->enc3,     sizeof(short),         snap->lin_size,   fp);
  fread(snap->nuc_enc,  sizeof(short),         snap->lin_size,   fp);
  fread(snap->nuc_enc5, sizeof(short),         n + 1,            fp);
  fread(snap->nuc_enc3, sizeof(short),         n + 1,            fp);
  fread(snap->ptype,    sizeof(char),          snap->ptype_size, fp);
  fread(snap->hc_mx,    sizeof(unsigned char), snap->hc_size,    fp);
  fread(snap->sequence, 1,                     n + 1,            fp);
  fread(snap->nuc_string, 1,                   n + 1,            fp);

  vrna_topk_cache_t *topk_cache = NULL;
  if (m == MODE_TOPK && K > 0) {
    size_t matrix_size = (size_t)(n + 1) * (size_t)(n + 2) / 2 + 2;

    topk_cache = (vrna_topk_cache_t *)calloc(1, sizeof(vrna_topk_cache_t));
    topk_cache->K = K;
    topk_cache->n = n;

    topk_cache->fml_splits = (unsigned short *)malloc(sizeof(unsigned short) * matrix_size * K);
    topk_cache->mb_splits  = (unsigned short *)malloc(sizeof(unsigned short) * matrix_size * K);
    topk_cache->int_k      = (unsigned short *)malloc(sizeof(unsigned short) * matrix_size * K);
    topk_cache->int_l      = (unsigned short *)malloc(sizeof(unsigned short) * matrix_size * K);
    topk_cache->fml_count  = (unsigned char *)malloc(sizeof(unsigned char) * matrix_size);
    topk_cache->mb_count   = (unsigned char *)malloc(sizeof(unsigned char) * matrix_size);
    topk_cache->int_count  = (unsigned char *)malloc(sizeof(unsigned char) * matrix_size);

    fread(topk_cache->fml_splits, sizeof(unsigned short), matrix_size * K, fp);
    fread(topk_cache->mb_splits,  sizeof(unsigned short), matrix_size * K, fp);
    fread(topk_cache->int_k,      sizeof(unsigned short), matrix_size * K, fp);
    fread(topk_cache->int_l,      sizeof(unsigned short), matrix_size * K, fp);
    fread(topk_cache->fml_count,  sizeof(unsigned char),  matrix_size,     fp);
    fread(topk_cache->mb_count,   sizeof(unsigned char),  matrix_size,     fp);
    fread(topk_cache->int_count,  sizeof(unsigned char),  matrix_size,     fp);
  }

  fclose(fp);

  cache_data_t *cd = (cache_data_t *)calloc(1, sizeof(cache_data_t));
  cd->mode       = m;
  cd->n          = n;
  cd->K          = K;
  cd->init_mfe   = (float)mfe_int / 100.0f;
  cd->sequence   = sequence;
  cd->structure  = structure;
  cd->snap       = snap;
  cd->topk_cache = topk_cache;

  return cd;
}


static void
cache_data_free(cache_data_t *cd)
{
  if (!cd) return;
  free(cd->sequence);
  free(cd->structure);
  snapshot_free(cd->snap);
  if (cd->topk_cache)
    vrna_topk_cache_free(cd->topk_cache);
  free(cd);
}


/* ================================================================
 * FASTA reader
 * ================================================================ */

static char *
read_sequence(const char *input)
{
  /* Check if input is a file */
  FILE *fp = fopen(input, "r");
  if (fp) {
    char line[MAX_SEQ];
    char *seq = (char *)calloc(MAX_SEQ, 1);
    size_t seq_len = 0;

    while (fgets(line, sizeof(line), fp)) {
      size_t len = strlen(line);
      while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
        line[--len] = '\0';

      if (line[0] == '>' || line[0] == ';')
        continue;

      if (len > 0) {
        memcpy(seq + seq_len, line, len);
        seq_len += len;
      }
    }
    seq[seq_len] = '\0';
    fclose(fp);

    /* Validate: must be RNA/DNA */
    for (size_t i = 0; i < seq_len; i++)
      seq[i] = toupper(seq[i]);

    return seq;
  }

  /* Not a file — treat as raw sequence */
  size_t len = strlen(input);
  char *seq = (char *)malloc(len + 1);
  for (size_t i = 0; i < len; i++)
    seq[i] = toupper(input[i]);
  seq[len] = '\0';
  return seq;
}


/* ================================================================
 * Mutation list parsing
 * ================================================================ */

typedef struct {
  unsigned int pos;
  char         base;
} mutation_t;


static int
parse_mutations_file(const char *path, mutation_t **out_mutations, int max_mutations)
{
  FILE *fp = fopen(path, "r");
  if (!fp) { fprintf(stderr, "Error: cannot open mutation file %s\n", path); return 0; }

  mutation_t *muts = (mutation_t *)malloc(sizeof(mutation_t) * max_mutations);
  int count = 0;
  char line[1024];

  while (fgets(line, sizeof(line), fp) && count < max_mutations) {
    /* Skip comments and empty lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
      continue;

    unsigned int pos;
    char base;
    if (sscanf(line, "%u %c", &pos, &base) == 2) {
      muts[count].pos  = pos;
      muts[count].base = toupper(base);
      count++;
    }
  }

  fclose(fp);
  *out_mutations = muts;
  return count;
}


/* ================================================================
 * Parallel fold worker
 * ================================================================ */

typedef struct {
  /* Input */
  cache_data_t  *cd;
  mutation_t    *mutations;
  int           n_mutations;
  int           worker_id;
  int           n_workers;

  /* Output */
  char          **results;
  int           n_results;
} worker_arg_t;


static void *
fold_worker(void *arg)
{
  worker_arg_t *wa = (worker_arg_t *)arg;
  cache_data_t *cd = wa->cd;

  vrna_fold_compound_t *fc = vrna_fold_compound(cd->sequence, NULL, VRNA_OPTION_MFE);
  char *work_struct = (char *)vrna_alloc(cd->n + 1);
  vrna_mfe(fc, work_struct);

  snapshot_restore(cd->snap, fc);

  /* Fix up the topk_cache indx pointer after deserialization */
  vrna_topk_cache_t *cache = cd->topk_cache;
  if (cache)
    cache->indx = fc->jindx;

  unsigned int flags = 0;
  if (cd->mode == MODE_TOPK)
    flags = VRNA_MFE_INC_TOPK;

  /* Allocate result buffer */
  int my_count = 0;
  for (int i = wa->worker_id; i < wa->n_mutations; i += wa->n_workers)
    my_count++;

  wa->results = (char **)malloc(sizeof(char *) * my_count);
  wa->n_results = 0;

  char *out_struct = (char *)vrna_alloc(cd->n + 1);

  for (int i = wa->worker_id; i < wa->n_mutations; i += wa->n_workers) {
    mutation_t *m = &wa->mutations[i];

    snapshot_restore(cd->snap, fc);

    double t0 = now_sec();
    float mfe = vrna_mfe_incremental_ex(fc, m->pos, m->base, flags, out_struct, cache);
    double elapsed = (now_sec() - t0) * 1e3;

    char orig = cd->sequence[m->pos - 1];
    const char *method = (cd->mode == MODE_TOPK) ? "topk" : "incremental";

    char *buf = (char *)malloc(cd->n + 256);
    snprintf(buf, cd->n + 256, "%u,%c,%c,%s,%.2f,%s,%.3f",
             m->pos, orig, m->base, method, mfe, out_struct, elapsed);

    wa->results[wa->n_results++] = buf;
  }

  free(out_struct);
  free(work_struct);
  vrna_fold_compound_free(fc);

  return NULL;
}


/* ================================================================
 * Subcommand: cache
 * ================================================================ */

static int
cmd_cache(int argc, char *argv[])
{
  const char *input  = NULL;
  const char *output = NULL;
  int K = 5;
  int mode = MODE_TOPK;
  int opt;

  /* Pre-scan for -inc/--inc (long-form) before getopt processes it */
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-inc") == 0 || strcmp(argv[i], "--inc") == 0) {
      mode = MODE_INC;
      argv[i] = (char *)"-i";  /* replace so getopt sees -i */
    }
  }

  optind = 1;
  while ((opt = getopt(argc, argv, "s:k:o:ih")) != -1) {
    switch (opt) {
      case 's': input = optarg; break;
      case 'k':
        K = atoi(optarg);
        mode = MODE_TOPK;
        break;
      case 'i': mode = MODE_INC; break;
      case 'o': output = optarg; break;
      case 'h':
        fprintf(stderr,
          "Usage: rapidrefold cache -s <input> [-k K | -inc] -o <output>\n"
          "\n"
          "  -s FILE/SEQ   Input FASTA file or raw RNA sequence\n"
          "  -k K          Build top-K cache with K candidates (default: 5)\n"
          "  -i / -inc     Exact incremental mode (no top-K, snapshot only)\n"
          "  -o FILE       Output binary cache file\n"
          "  -h            Show this help\n");
        return 0;
      default: return 1;
    }
  }

  if (!input) { fprintf(stderr, "Error: -s (input) required\n"); return 1; }
  if (!output) { fprintf(stderr, "Error: -o (output) required\n"); return 1; }

  char *seq = read_sequence(input);
  if (!seq || strlen(seq) == 0) {
    fprintf(stderr, "Error: could not read sequence from '%s'\n", input);
    free(seq); return 1;
  }

  unsigned int n = strlen(seq);
  fprintf(stderr, "Sequence length: %u\n", n);
  fprintf(stderr, "Mode: %s", mode == MODE_TOPK ? "top-K" : "incremental");
  if (mode == MODE_TOPK) fprintf(stderr, " (K=%d)", K);
  fprintf(stderr, "\n");

  /* Initial fold */
  fprintf(stderr, "Folding...\n");
  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE);
  char *structure = (char *)vrna_alloc(n + 1);

  double t0 = now_sec();
  float init_mfe = vrna_mfe(fc, structure);
  double t_fold = (now_sec() - t0) * 1e3;

  fprintf(stderr, "  MFE: %.2f kcal/mol (%.1f ms)\n", init_mfe, t_fold);

  /* Save DP snapshot */
  dp_snapshot_t *snap = snapshot_alloc(n);
  snapshot_save(snap, fc);

  /* Build top-K cache if needed */
  vrna_topk_cache_t *topk_cache = NULL;
  if (mode == MODE_TOPK) {
    fprintf(stderr, "Building top-K cache (K=%d)...\n", K);
    t0 = now_sec();
    topk_cache = vrna_topk_cache_from_fold(fc, K);
    double t_cache = (now_sec() - t0) * 1e3;
    fprintf(stderr, "  Cache built in %.1f ms\n", t_cache);
  }

  /* Write binary cache */
  fprintf(stderr, "Writing cache to %s...\n", output);
  int rc = write_cache_file(output, mode, seq, structure, init_mfe, K, snap, topk_cache);

  if (rc == 0)
    fprintf(stderr, "Done.\n");

  /* Cleanup */
  if (topk_cache) vrna_topk_cache_free(topk_cache);
  snapshot_free(snap);
  vrna_fold_compound_free(fc);
  free(structure);
  free(seq);

  return rc;
}


/* ================================================================
 * Subcommand: fold
 * ================================================================ */

static int
cmd_fold(int argc, char *argv[])
{
  const char *cache_path = NULL;
  const char *output     = NULL;
  const char *mut_file   = NULL;
  int n_workers = 1;

  mutation_t *cmd_mutations = NULL;
  int cmd_n_mutations = 0;
  int cmd_mut_cap = 64;

  cmd_mutations = (mutation_t *)malloc(sizeof(mutation_t) * cmd_mut_cap);

  /* Manual argument parsing to handle: fold <cache> -m POS BASE ...
   * argv[0] is the subcommand name "fold", skip it */
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      fprintf(stderr,
        "Usage: rapidrefold fold <cache.bin> [options]\n"
        "\n"
        "  <cache.bin>    Binary cache file from 'rapidrefold cache'\n"
        "  -m POS BASE    Single mutation (1-based position, nucleotide)\n"
        "                 Can be specified multiple times\n"
        "  -M FILE        Mutation list file (one 'POS BASE' per line)\n"
        "  -j N           Number of parallel workers (default: 1)\n"
        "  -o FILE        Output CSV file (default: stdout)\n"
        "  -h             Show this help\n"
        "\n"
        "Mutation file format:\n"
        "  # Lines starting with '#' are comments\n"
        "  5 C\n"
        "  10 U\n"
        "  23 A\n");
      free(cmd_mutations);
      return 0;
    } else if (strcmp(argv[i], "-m") == 0 && i + 2 < argc) {
      unsigned int pos = atoi(argv[i + 1]);
      char base = toupper(argv[i + 2][0]);
      if (cmd_n_mutations >= cmd_mut_cap) {
        cmd_mut_cap *= 2;
        cmd_mutations = (mutation_t *)realloc(cmd_mutations, sizeof(mutation_t) * cmd_mut_cap);
      }
      cmd_mutations[cmd_n_mutations].pos = pos;
      cmd_mutations[cmd_n_mutations].base = base;
      cmd_n_mutations++;
      i += 2;
    } else if (strcmp(argv[i], "-M") == 0 && i + 1 < argc) {
      mut_file = argv[++i];
    } else if (strcmp(argv[i], "-j") == 0 && i + 1 < argc) {
      n_workers = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      output = argv[++i];
    } else if (argv[i][0] != '-' && !cache_path) {
      cache_path = argv[i];
    }
  }

  if (!cache_path) {
    fprintf(stderr, "Error: cache file required\n");
    free(cmd_mutations);
    return 1;
  }

  /* Load mutations from file if specified */
  mutation_t *file_mutations = NULL;
  int file_n_mutations = 0;

  if (mut_file) {
    file_n_mutations = parse_mutations_file(mut_file, &file_mutations, 100000);
    if (file_n_mutations == 0) {
      fprintf(stderr, "Warning: no mutations read from %s\n", mut_file);
    }
  }

  /* Merge mutations */
  int total_mutations = cmd_n_mutations + file_n_mutations;
  mutation_t *all_mutations = (mutation_t *)malloc(sizeof(mutation_t) * (total_mutations + 1));
  if (cmd_n_mutations > 0)
    memcpy(all_mutations, cmd_mutations, sizeof(mutation_t) * cmd_n_mutations);
  if (file_n_mutations > 0)
    memcpy(all_mutations + cmd_n_mutations, file_mutations, sizeof(mutation_t) * file_n_mutations);

  free(cmd_mutations);
  free(file_mutations);

  if (total_mutations == 0) {
    fprintf(stderr, "Error: no mutations specified (-m or -M required)\n");
    free(all_mutations);
    return 1;
  }

  /* Load cache */
  fprintf(stderr, "Loading cache from %s...\n", cache_path);
  cache_data_t *cd = read_cache_file(cache_path);
  if (!cd) {
    free(all_mutations);
    return 1;
  }

  fprintf(stderr, "  Sequence: %u nt, mode: %s",
          cd->n, cd->mode == MODE_TOPK ? "top-K" : "incremental");
  if (cd->mode == MODE_TOPK)
    fprintf(stderr, " (K=%u)", cd->K);
  fprintf(stderr, "\n");
  fprintf(stderr, "  Mutations to process: %d\n", total_mutations);

  /* Validate mutations */
  for (int i = 0; i < total_mutations; i++) {
    if (all_mutations[i].pos < 1 || all_mutations[i].pos > cd->n) {
      fprintf(stderr, "Error: mutation position %u out of range [1, %u]\n",
              all_mutations[i].pos, cd->n);
      cache_data_free(cd);
      free(all_mutations);
      return 1;
    }
    char b = toupper(all_mutations[i].base);
    if (b != 'A' && b != 'C' && b != 'G' && b != 'U' && b != 'T') {
      fprintf(stderr, "Error: invalid base '%c' at mutation %d\n",
              all_mutations[i].base, i + 1);
      cache_data_free(cd);
      free(all_mutations);
      return 1;
    }
    all_mutations[i].base = (b == 'T') ? 'U' : b;
  }

  /* Clamp workers */
  if (n_workers < 1) n_workers = 1;
  if (n_workers > total_mutations) n_workers = total_mutations;

  /* Run workers */
  fprintf(stderr, "Folding with %d worker(s)...\n", n_workers);
  double t_start = now_sec();

  worker_arg_t *workers = (worker_arg_t *)calloc(n_workers, sizeof(worker_arg_t));
  pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * n_workers);

  for (int w = 0; w < n_workers; w++) {
    workers[w].cd = cd;
    workers[w].mutations = all_mutations;
    workers[w].n_mutations = total_mutations;
    workers[w].worker_id = w;
    workers[w].n_workers = n_workers;
  }

  if (n_workers == 1) {
    fold_worker(&workers[0]);
  } else {
    for (int w = 0; w < n_workers; w++)
      pthread_create(&threads[w], NULL, fold_worker, &workers[w]);
    for (int w = 0; w < n_workers; w++)
      pthread_join(threads[w], NULL);
  }

  double t_total = (now_sec() - t_start) * 1e3;

  /* Write output */
  FILE *out = output ? fopen(output, "w") : stdout;
  if (!out) {
    fprintf(stderr, "Error: cannot open %s for writing\n", output);
    /* cleanup below */
  } else {
    fprintf(out, "pos,orig,new,method,mfe,structure,time_ms\n");

    /* Collect results in order */
    for (int i = 0; i < total_mutations; i++) {
      int w = i % n_workers;
      int local_idx = i / n_workers;
      if (local_idx < workers[w].n_results)
        fprintf(out, "%s\n", workers[w].results[local_idx]);
    }

    if (output) fclose(out);
  }

  fprintf(stderr, "  Total: %.1f ms (%.2f ms/mutation)\n",
          t_total, t_total / total_mutations);

  /* Cleanup */
  for (int w = 0; w < n_workers; w++) {
    for (int r = 0; r < workers[w].n_results; r++)
      free(workers[w].results[r]);
    free(workers[w].results);
  }
  free(workers);
  free(threads);
  free(all_mutations);
  cache_data_free(cd);

  return 0;
}


/* ================================================================
 * Main
 * ================================================================ */

static void
print_usage(void)
{
  fprintf(stderr,
    "RapidRefold — Fast incremental RNA refolding under point mutations\n"
    "\n"
    "Usage:\n"
    "  rapidrefold cache -s <input> [-k K | -inc] -o <output>\n"
    "  rapidrefold fold  <cache.bin> [-m POS BASE | -M file] [-j N] [-o out.csv]\n"
    "\n"
    "Subcommands:\n"
    "  cache    Fold a sequence and save DP state to a binary cache file\n"
    "  fold     Load a cache and refold under specified mutations\n"
    "\n"
    "Run 'rapidrefold <subcommand> -h' for detailed help.\n");
}


int
main(int argc, char *argv[])
{
  if (argc < 2) {
    print_usage();
    return 1;
  }

  if (strcmp(argv[1], "cache") == 0) {
    return cmd_cache(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "fold") == 0) {
    return cmd_fold(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    print_usage();
    return 0;
  } else {
    fprintf(stderr, "Error: unknown subcommand '%s'\n\n", argv[1]);
    print_usage();
    return 1;
  }
}
