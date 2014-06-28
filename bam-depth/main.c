#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "bedidx.h"

typedef BGZF bamFile;
typedef hts_itr_t *bam_iter_t;
typedef bam_hdr_t bam_header_t;
typedef hts_idx_t bam_index_t;

typedef struct {     // auxiliary data structure
    bamFile *fp;      // the file handler
    bam_iter_t iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;


// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret = aux->iter? hts_itr_next(aux->fp, aux->iter, b, (hts_readrec_f)(bam_readrec), 0) : bam_read1(aux->fp, b);
    if (!(b->core.flag&BAM_FUNMAP)) {
        if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
        else if (aux->min_len && bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
    }
    return ret;
}

int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *beg, int *end)
{
    const char *name_lim = hts_parse_reg(str, beg, end);
    char *name = malloc(name_lim - str + 1);
    memcpy(name, str, name_lim - str);
    name[name_lim - str] = '\0';
    *ref_id = bam_name2id(header, name);
    free(name);
    if (*ref_id == -1) return -1;
    return *beg <= *end? 0 : -1;
}

int main(int argc, char *argv[])
{
    int i, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, avecov = 0, window_step_size = 0, window_size = 0;
    const bam_pileup1_t **plp;
    //int **window_matrix;   // matrix for the coverage in the window for each of the bam files
    char *reg = 0; // specified region
    void *bed = 0; // BED data structure
    bam_header_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;

    // parse the command line
    while ((n = getopt(argc, argv, "r:b:q:Q:l:f:aw:s:")) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
            case 'b': bed = bed_read(optarg); break; // BED or position list file can be parsed now
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'a': avecov = 1; break;
            case 'w': window_size = atoi(optarg); break;
            //case 's': window_step_size = atoi(optarg); break;
        }
    }
    if (optind == argc || (avecov && bed )) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: samtools depth [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -a                  Output average coverage for regions (default: per base)\n");
        fprintf(stderr, "   -b <bed>            list of positions or regions (not compatible with -a)\n");
        fprintf(stderr, "   -l <int>            minQLen\n");
        fprintf(stderr, "   -q <int>            base quality threshold\n");
        fprintf(stderr, "   -Q <int>            mapping quality threshold\n");
        fprintf(stderr, "   -r <chr:from-to>    region\n");
        fprintf(stderr, "   -w <int>            window size\n");
        //fprintf(stderr, "   -s <int>            window step size\n");
        fprintf(stderr, "\n");
        return 1;
    }

    // initialize the auxiliary data structures
    n = argc - optind; // the number of BAMs on the command line
    data = calloc(n, sizeof(void*)); // data[i] for the i-th input
    beg = 0; end = 1<<30; tid = -1;  // set the default region
    for (i = 0; i < n; ++i) {
        bam_header_t *htmp;
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(argv[optind+i], "r"); // open BAM
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        htmp = bam_hdr_read(data[i]->fp);         // read the BAM header
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
            if (reg) bam_parse_region(h, reg, &tid, &beg, &end); // also parse the region
        } else bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header
        if (tid >= 0) { // if a region is specified and parsed successfully
            bam_index_t *idx = bam_index_load(argv[optind+i]);  // load the index
            data[i]->iter = bam_itr_queryi(idx, tid, beg, end); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
        }
    }
    /*if(window_size) {
        if(window_step_size == 0) {
            window_step_size = window_size;
        }
        window_matrix = calloc(n, sizeof(int*));
        for(int i = 0; i < n; ++i) {
            window_matrix[i] = calloc(window_size, sizeof(int));
            for(int j = 0; j < window_size; ++j) {
                window_matrix[i][j] = 0;
            }
        }
    }*/

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
    int *total_reads = calloc(n, sizeof(int));  // store the total number of reads that mapped per chrom per bam
    // initialize
    for (i = 0; i < n; ++i) {
        total_reads[i] = 0;
    }
    int prev_tid = -1;  // the id of the previous positions tid
    int windows_in_chrom = 1;
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip

        if(avecov || window_size) {
            if(prev_tid != -1 && tid != prev_tid) {  // we've arrived at a new chrom
                int start = beg;
                float reg_length;
                if(reg) {
                    reg_length = (float)(end+1 - beg);
                } else if (window_size) {
                    reg_length = (float)(h->target_len[prev_tid] % window_size);
                    start = h->target_len[prev_tid] - (int)reg_length;
                } else {
                    reg_length = (float)h->target_len[prev_tid];
                }
                printf("%s\t%d\t%.0f", h->target_name[prev_tid], start, reg_length);
                for (i = 0; i < n; ++i) {
                    printf("\t%.2f", (float)total_reads[i]/reg_length);
                    total_reads[i] = 0;  // reset for next chrom
                }
                putchar('\n');
            }
            if(window_size && pos > window_size * windows_in_chrom) {
                printf("%s\t%d\t%d", h->target_name[prev_tid],0, windows_in_chrom*window_size);
                for (i = 0; i < n; ++i) {
                    printf("\t%.2f", (float)total_reads[i]/(float)window_size);
                    total_reads[i] = 0;  // reset for next window
                }
                putchar('\n');
                ++windows_in_chrom;
            }
        }else {
            printf("%s\t%d", h->target_name[tid], pos+1);
        }
        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j, m = 0;
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
                else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
            }
            if(avecov || window_size) {
                total_reads[i] += n_plp[i] - m;  // add this positions depth to the total
            } else {
                printf("\t%d", n_plp[i] - m); // this the depth to output
            }
        }
        if(!avecov) {
            putchar('\n');
        }
        prev_tid = tid;
    }
    if(avecov) {
        if(prev_tid != -1) {  // last chrom
            fputs(h->target_name[prev_tid], stdout);
            float reg_length;
            if(reg) {
                reg_length = (float)(end+1 - beg);
            } else {
                reg_length = (float)h->target_len[prev_tid];
            }
            printf("\t%.0f", reg_length);
            for (i = 0; i < n; ++i) {
                printf("\t%.2f", (float)total_reads[i]/reg_length);
            }
            putchar('\n');
        }
    }
    free(n_plp); free(plp); free(total_reads);
    bam_mplp_destroy(mplp);

    bam_hdr_destroy(h);
    for (i = 0; i < n; ++i) {
        bgzf_close(data[i]->fp);
        if (data[i]->iter) bam_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(reg);
    if (bed) bed_destroy(bed);

    /*if(window_size) {
        for(int i = 0; i < n; ++i) {
            free(window_matrix[i]);
        }
        free(window_matrix);
    }*/
    return 0;
}
