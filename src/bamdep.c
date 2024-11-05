#include "bamdep.h"

int main(int argc, char *argv[])
{
	if (argc == 1)
		usage();
	setenv("FONTCONFIG_PATH", "/etc/fonts", 1);
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->mpq = MINMQ;
	arg->len = MINQL;
	prs_arg(argc, argv, arg);
	samFile *fp = sam_open(arg->in, "r");
	if (!fp)
		error("Error: failed to read input bam [%s]\n", arg->in);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	int i, ci = -1;
	// check ctg
	if (arg->ctg)
	{
		for (i = 0; i < hdr->n_targets; ++i)
		{
			if (!strcmp(hdr->target_name[i], arg->ctg))
			{
				ci = i;
				break;
			}
		}
		if (ci == -1)
		{
			fprintf(stderr, "%s Error: specified contig [%s] not found in bam\n", ERR, arg->ctg);
			fprintf(stderr, "Please use samtools view -H %s to view contigs contained\n", arg->in);
			exit(EXIT_FAILURE);
		}
	}
	// multiple ctg offsets
	uint64_t gl = 0, nd = 0; // genome length in total
	kh_t *os = kh_init();
	ld_os(hdr, ci, os, &gl);
	/* dbg os
	for (i = 0; i < hdr->n_targets; ++i)
		printf("%s\t%d\t%"PRIu64"\n", hdr->target_name[i], hdr->target_len[i], kh_xval(os, i));
	printf("%"PRIu64"\n", gl);
	*/
	uint32_t md = 0;
	dp_t *dp = NULL;
	ld_dp(arg->in, arg->ctg, &dp, &md, &nd);
	/* debug dp
	for (i = 0; i < nd; ++i)
		printf("%s\t%d\t%d\t%d\n", hdr->target_name[dp[i].tid], dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
	*/
	// dep plot
	cairo_t *cr = NULL;
	cairo_surface_t *sf = NULL;
	if (ends_with(arg->out, ".svg"))
		sf = cairo_svg_surface_create(arg->out, WIDTH * 1.02, HEIGHT);
	else if (ends_with(arg->out, ".png"))
		sf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, WIDTH, HEIGHT);
	else
		error("Error: unsupported plot format: [%s]\n", arg->out);
	cr = cairo_create(sf);
	char *tt, an[NAME_MAX] = {'\0'};
	// prepare title
	char *p = strrchr(arg->in, '/');
	asprintf(&tt, "%s", p ? p + 1 : arg->in);
	*strstr(tt, ".bam") = '\0';
	if (ci != -1)
		asprintf(&tt, "%s: %s", tt, arg->ctg);
	// prepare anno
	prep_an(dp, nd, gl, an);
	draw_canvas(sf, cr, hdr, ci, os, tt, arg->sub, an, md, gl);
	cairo_save(cr);
	cairo_scale(cr, DIM_X, DIM_Y);
	// iterate dp array and draw to cr
	for (i = 0; i < nd; ++i)
		draw_ped1(cr, os, md, gl, dp + i);
	cairo_restore(cr);
	draw_axis(cr, md, gl);
	free(dp);
	free(tt);
	bam_hdr_destroy(hdr);
	hts_close(fp);
	if (ends_with(arg->out, ".png"))
		cairo_surface_write_to_png(cairo_get_target(cr), arg->out);
	cairo_surface_destroy(sf);
	cairo_destroy(cr);
	kh_destroy(os);
	free(arg);
	return 0;
}

void ld_os(bam_hdr_t *hdr, int ci, kh_t *os, uint64_t *gl)
{
	int i;
	uint64_t shift = 0;
	if (ci != -1)
	{
		kh_ins(os, ci, 0);
		*gl = hdr->target_len[ci];
	}
	else
	{
		for (i = 0; i < hdr->n_targets; ++i)
		{
			kh_ins(os, i, shift);
			shift += hdr->target_len[i];
		}
		*gl = shift;
	}
}

void prep_an(const dp_t *dp, uint64_t nd, uint64_t gl, char *an)
{
	uint64_t i;
	double cl = 0.0f;
	for (i = 0; i < nd; ++i)
		cl += dp[i].len;
	cl /= gl / 100;
	if (cl == 0)
		snprintf(an, NAME_MAX, "Genome coverage: N/A");
	else if (cl < 0.001)
		snprintf(an, NAME_MAX, "Genome coverage: <0.01‰");
	else if (cl < 0.01)
		snprintf(an, NAME_MAX, "Genome coverage: <0.01%%");
	else
		snprintf(an, NAME_MAX, "Genome coverage: %.*f%%", fmin(100, cl) == 100 ? 0 : 2, fmin(100, cl));
}

void draw_canvas(cairo_surface_t *sf, cairo_t *cr, bam_hdr_t *hdr, int ci,
		const kh_t *os, const char *tt, const char *st, const char *an, uint32_t md,
		uint64_t gl)
{
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
	draw_rrect(cr);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_translate(cr, MARGIN / 1.25, MARGIN / 2.0);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	if (tt && strlen(tt)) // title
	{
		cairo_set_font_size(cr, 24.0);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		cairo_text_extents(cr, tt, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = (st && strlen(st)) ? -ext.height / 2 + ext.y_bearing * 1.25 : ext.height + ext.y_bearing * 2.5;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, tt);
	}
	if (st && strlen(st)) // sub-title
	{
		cairo_set_font_size(cr, 20.0);
		cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_ITALIC, CAIRO_FONT_WEIGHT_BOLD);
		cairo_text_extents(cr, st, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.height + ext.y_bearing * 2;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, st);
	}
	if (an && strlen(an)) // zlab
	{
		char zlab[NAME_MAX];
		sprintf(zlab, "%s", an);
		cairo_set_font_size(cr, 18.0);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, zlab, &ext);
		x = DIM_X - (ext.width + ext.x_bearing);
		y = ext.height / 2 + ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, zlab);
	}
	// xlab
	char xlab[] = "Genome coordinates";
	char ylab[] = "Depth (PE)";
	cairo_set_font_size(cr, 18.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, xlab);
	cairo_save(cr);
	// ylab
	cairo_translate(cr, MARGIN / 1.25, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, ylab, &ext);
	cairo_move_to(cr, MARGIN / 2.5, -MARGIN * 1.25);
	cairo_show_text(cr, ylab);
	cairo_restore(cr);
	// contig shades
	int i;
	cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
	if (ci == -1)
	{
		for (i = 1; i < hdr->n_targets; i += 2)
		{
			x = kh_xval(os, i);
			cairo_rectangle(cr, (double)DIM_X * x / gl, 0, (double)DIM_X * hdr->target_len[i] / gl, DIM_Y);
			cairo_fill(cr);
		}
	}
}

void draw_axis(cairo_t *cr, uint32_t md, uint64_t gl)
{
	double x, y;
	cairo_text_extents_t ext;
	draw_arrow(cr, 0, DIM_Y, DIM_X, DIM_Y); // xaxis
	double w1 = 1.0, w2 = 1.0;
	cairo_set_font_size(cr, 16.0);
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, 0, 0);
	cairo_line_to(cr, 0, DIM_Y); // yaxis
	draw_yticks(cr, md);
	char buf[sizeof(uint64_t) * 8 + 1];
	if (gl >= 1e9)
		sprintf(buf, "%.2fG", gl * 1.0e-9);
	else if (gl >= 1e6)
		sprintf(buf, "%.2fM", gl * 1.0e-6);
	else if (gl >= 1e3)
		sprintf(buf, "%.2fK", gl * 1.0e-3);
	else
		sprintf(buf, "%"PRIu64, gl);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, buf, &ext);
	x = DIM_X - ext.width - ext.x_bearing;
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, buf);
}

void draw_ped1(cairo_t *cr, kh_t *os, uint32_t md, uint64_t gl, dp_t *dp)
{
	double w1 = 1.0, w2 = 1.0, x, y, w, h;
	cairo_device_to_user_distance(cr, &w1, &w2);
	double lw = fmin(w1, w2) / 2;
	cairo_set_line_width(cr, lw);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	double ymx = ceil(log10(md)) + 1;
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	if (dp->len <= 5) // use hist instead of rectangle to make it visible
	{
		cairo_set_line_width(cr, lw);
		x = (double)(dp->pos + kh_xval(os, dp->tid) + (double)dp->len / 2) / gl;
		y = 1 - (log10(dp->dep) + 1) / ymx;
		cairo_move_to(cr, x, 1);
		cairo_line_to(cr, x, y);
		cairo_stroke(cr);
	}
	else
	{
		x = (double)(dp->pos + kh_xval(os, dp->tid)) / gl;
		y = 1 - (log10(dp->dep) + 1) / ymx;
		w = (double)dp->len / gl;
		h = (log10(dp->dep) + 1) / ymx;
		cairo_rectangle(cr, x, y, w, h);
		cairo_fill(cr);
	}
}

int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	int ret;
	while (true)
	{
		ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if (ret < 0)
			break;
		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if ((int)b->core.qual < aux->min_mapQ)
			continue;
		if (aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len)
			continue;
		break;
	}
	return ret;
}

void ld_dp(const char *fn, const char *ctg, dp_t **dp, uint32_t *md, uint64_t *nd)
{
	int64_t d = 0, m = CHUNK, n, p = 0;
	*dp = malloc(m * sizeof(dp_t));
	int tid, pos, beg = 0, end = INT_MAX, n_plp;
	aux_t *data = calloc(1, sizeof(aux_t));
	data->fp = hts_open(fn, "r");
	data->hdr = sam_hdr_read(data->fp);
	hts_idx_t *idx = sam_index_load(data->fp, fn);
	if (ctg)
		data->iter = sam_itr_querys(idx, data->hdr, ctg);
	hts_idx_destroy(idx);
	bam_hdr_t *h = data->hdr;
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**)&data);
	bam_mplp_set_maxcnt(mplp, MAXDP);
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t));
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0)
	{
		if (pos < beg)
			continue;
		if (pos >= end || tid >= h->n_targets)
			break;
		if ((n = n_plp))
		{
			if (n != d)
			{
				if (p)
				{
					(*dp)[*nd].len = pos - (*dp)[*nd].pos + 1;
					d = p = 0;
					if (++*nd == m)
					{
						m <<= 1;
						*dp = realloc(*dp, m * sizeof(dp_t));
					}
				}
				(*dp)[*nd].tid = tid;
				(*dp)[*nd].pos = pos;
				(*dp)[*nd].dep = n;
				*md = fmax(*md, n);
			}
			d = n;
			p = pos;
		}
		else
		{
			if (d)
			{
				(*dp)[*nd].len = pos - (*dp)[*nd].pos + 1;
				if (++*nd == m)
				{
					m <<= 1;
					*dp = realloc(*dp, m * sizeof(dp_t));
				}
			}
			d = p = 0;
		}
	}
	if (d)
	{
		(*dp)[*nd].len = pos - (*dp)[*nd].pos + 1;
		if (++*nd == m)
		{
			m <<= 1;
			*dp = realloc(*dp, m * sizeof(dp_t));
		}
	}
	d = p = 0;
	*dp = realloc(*dp, *nd * sizeof(dp_t));
	bam_hdr_destroy(data->hdr);
	if (data->fp)
		hts_close(data->fp);
	hts_itr_destroy(data->iter);
	free(data);
}

void dump_dp(bam_hdr_t *hdr, dp_t *dp, uint64_t nd, const char *out)
{
	uint64_t i;
	kstring_t ks = {0, 0, NULL};
	if (out && ends_with(out, ".gz"))
	{
		BGZF *fp = bgzf_open(out, "w");
		for (i = 0; i < nd; ++i)
		{
			ksprintf(&ks, "%s\t%"PRIu64"\t%"PRIu64"\t%d\n", hdr->target_name[dp[i].tid],
					dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
			if (ks.l >= BGZF_BLOCK_SIZE)
			{
				if (ks.l != bgzf_write(fp, ks.s, ks.l))
					error("Error writing to file [%s]\n", out);
				ks.l = 0;
			}
		}
		if (ks.l && ks.l != bgzf_write(fp, ks.s, ks.l))
			error("Error writing to file [%s]\n", out);
		ks.l = 0;
		bgzf_close(fp);
		tbx_conf_t conf = tbx_conf_bed;
		if(tbx_index_build3(out, NULL, 14, 8, &conf))
			error("Error building index for [%s]\n", out);
	}
	else
	{
		FILE *fp = out ? fopen(out, "w") : stdout;
		for (i = 0; i < nd; ++i)
			fprintf(fp, "%s\t%"PRIu64"\t%"PRIu64"\t%d\n", hdr->target_name[dp[i].tid],
					dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
		fclose(fp);
	}
	ks_release(&ks);
}

int is_gzip(const char *fn)
{
	char buf[2];
	FILE *fp;
	int gzip = 0;
	if ((fp = fopen(fn, "rb")) == NULL)
		error("[ERROR] Unable to open file: %s\n", fn);
	if (fread(buf, 1, 2, fp) == 2)
		if (((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b))
			gzip = 1;
	fclose(fp);
	return gzip;
}

bool ends_with(const char *str, const char *sfx)
{
	int ret = 0;
	int str_len = strlen(str);
	int sfx_len = strlen(sfx);
	if ((str_len >= sfx_len) && (0 == strcasecmp(str + (str_len-sfx_len), sfx)))
		ret = 1;
	return ret;
}

void prs_arg(int argc, char **argv, arg_t *arg)
{
	int c = 0;
	ketopt_t opt = KETOPT_INIT;
	const char *opt_str = "i:o:m:l:s:c:hv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 'm': arg->mpq = atoi(opt.arg); break;
			case 'l': arg->len = atoi(opt.arg); break;
			case 's': arg->sub = opt.arg; break;
			case 'c': arg->ctg = opt.arg; break;
			case 'h': usage(); break;
			case 'v':
				if (strlen(BRANCH_COMMIT))
					printf("%s (%s)\n", VERSION, BRANCH_COMMIT);
				else
					puts(VERSION);
				exit(EXIT_SUCCESS);
			case '?':
				printf("Invalid option: [%c]\n", opt.opt); exit(EXIT_SUCCESS);
			default:
				printf("Invalid option: [%s]", opt.arg); exit(EXIT_SUCCESS);
		}
	}
	if (!arg->in || access(arg->in, R_OK))
		error("Error: input bam is unspecified or inaccessible!\n");
	if (!(ends_with(arg->in, ".bam") && is_gzip(arg->in)))
		error("Oops! only bam input is supported. Invalid bam file [%s]\n", arg->in);
	char bai[PATH_MAX];
	snprintf(bai, PATH_MAX, "%s.bai", arg->in);
	if (access(bai, R_OK))
		error("Error: bam's index file (.bai) is required, please use samtools sort and index to create it.\n");
}

int strlen_wo_esc(const char *str)
{
	int length = 0;
	while (*str != '\0')
	{
		// Check if this is the start of an ANSI escape sequence (ESC [ ... m)
		if (*str == '\033' && *(str + 1) == '[')
		{
			// Move past \033[
			str += 2;
			// Skip until we reach 'm', which ends the ANSI sequence
			while (*str != '\0' && *str != 'm')
				str++;
			// Move past 'm' if we found it
			if (*str == 'm')
				str++;
			continue;
		}
		// Count visible characters
		length++;
		str++;
	}
	return length;
}

void horiz(int n, bool bold)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	n = (w.ws_col >= n) ? n : w.ws_col;
	printf("\e[%s90m", bold ? "1;" : "");
	while (n--) fputs("\xe2\x94\x80", stdout);
	printf("\e[%s0m\n", bold ? "0;" : "");
}

void dep_sch()
{
	puts("\e[90m                            ┌──┐\n"
		"                    ┌─┐    ┌┘  └─┐\n"
		"              ──────┴─┴────┴─────┴───────\e[0m\n"
		"                    \e[31m——→\e[0m    \e[34m←——\e[0m\n"
		"                            \e[31m——\e[0m\e[34m←\e[0m\e[31m→\e[0m\e[34m——\e[0m");
}

void usage()
{
	int w = 58;
	horiz(w, true);
	char title[] = "\e[1mCalculate the depth of coverage of microbial genome\e[0m";
	int title_len = strlen_wo_esc(title);
	printf("%*.*s\n", (int)((w - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	dep_sch();
	horiz(w, false);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m --in <bam> --out <bed>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m   Input BAM file with bai index");
	puts("  -o, --out \e[3mSTR\e[0m    Output depth plot png \e[90m[${prefix}.png]\e[0m");
	puts("  -s, --sub \e[3mFILE\e[0m   Sub-title of depth plot \e[90m[none]\e[0m");
	puts("  -c, --ctg \e[3mSTR\e[0m    Restrict analysis to this contig \e[90m[none]\e[0m");
	putchar('\n');
	puts("  -h               Show help message");
	puts("  -v               Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz(w, true);
	exit(1);
}
