#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>
#include <sys/ioctl.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "ketopt.h"
#include "dplot.h"
#include "thpool.h"
#include "version.h"

extern const char *__progname;

#define VERSION "0.1.0"
#define CHUNK 0xFFFF
#define MAXDP 0xFFFFF
#define MINQL 35
#define MINMQ 10

#define ARR "\e[90m\xE2\x97\x82\e[0m"
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define ERR "\e[31m\xE2\x9C\x97\e[0m"
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

#define error(msg, ...) do {                                 \
    fputs("\e[31m\xE2\x9C\x97\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
    exit(EXIT_FAILURE);                                      \
} while (0)

// argument struct
typedef struct
{
	char *in, *out, *sub, *ann, *ctg, *dep, *dup;
	int mpq, len;
	bool bold;
} arg_t;

typedef struct
{
	char *in, *ctg;
	dp_t **dp;
	uint32_t *md;
	uint64_t *nd;
	bool dup;
} op_t;

typedef struct
{
	bam_hdr_t *hdr;
	dp_t *dp;
	uint64_t nd;
	char *out;
} dd_t;

static ko_longopt_t long_options[] = {
	{ "in",        ko_required_argument, 'i' },
	{ "out",       ko_required_argument, 'o' },
	{ "mpq",       ko_required_argument, 'm' },
	{ "len",       ko_required_argument, 'l' },
	{ "sub",       ko_required_argument, 's' },
	{ "ctg",       ko_required_argument, 'c' },
	{ "dep",       ko_required_argument, 'd' },
	{ "dup",       ko_required_argument, 'D' },
	{ "bold",      ko_no_argument, 'B' },
	{ "help",      ko_no_argument, 'h' },
	{ "version",   ko_no_argument, 'v' },
	{ NULL, 0, 0 }
};

typedef struct // auxiliary data structure
{
	samFile *fp;     // the file handle
	bam_hdr_t *hdr;  // the file header
	hts_itr_t *iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
	bool keep_dup; // keep duplicates
} aux_t;

void prs_arg(int argc, char **argv, arg_t *arg);
void ld_os(bam_hdr_t *hdr, int ci, kh_t *os, uint64_t *gl);
bool bam_has_dup(const char *fn);
int read_bam(void *data, bam1_t *b);
void ld_dp(void *op);
void dump_dp(void *op);
void prep_an(const dp_t *dp, uint64_t nd, uint64_t gl, char *an);
void draw_canvas(cairo_surface_t *sf, cairo_t *cr, bam_hdr_t *hdr, int ci,
		const kh_t *os, const char *tt, const char *st, const char *an, uint32_t md,
		uint64_t gl);
void draw_ped1(cairo_t *cr, kh_t *os, uint32_t md, uint64_t gl, bool bold, bool dup,
		bool svg, dp_t *dp);
void draw_axis(cairo_t *cr, uint32_t md, const char *ctg, uint32_t n_targets,
		uint64_t gl);
int is_gzip(const char *fn);
bool ends_with(const char *str, const char *sfx);
int strlen_wo_esc(const char *str);
void horiz(int n, bool bold);
void usage();
