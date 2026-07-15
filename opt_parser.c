#include "opt_parser.h"


// long-only option identifiers (kept out of the short-option ASCII range)
enum {
    OPT_BLEND = 1000,
    OPT_MINIMIZER,
    OPT_SYNCMER,
    OPT_LCP,
    OPT_DCT
};

int file_exists(char *filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

static const char *sketch_type_name(sketch_type_t t) {
    switch (t) {
        case SKETCH_BLEND:     return "blend";
        case SKETCH_MINIMIZER: return "minimizer";
        case SKETCH_SYNCMER:   return "syncmer";
        case SKETCH_LCP:       return "lcp";
        default:               return "unknown";
    }
}

void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Options:\n"
        "  -f, --fasta <file>        FASTA file\n"
        "  -q, --fastq <file>        FASTQ file\n"
        "  -o, --output <file>       BED file\n"
        "  -k, --kmer <int>          k-mer size [%d]\n"
        "  -w, --window <int>        window size [%d]\n"
        "  -b, --blend-bits <int>    number of bits for hash [%d]\n"
        "  -n, --n-neighbors <int>   number of neighbour [%d]\n"
        "  -s, --smer-size <int>     smmer size [%d]\n"
        "  -l, --lcp-level <int>     lcp level [%d]\n"
        "      --dct-count <int>     lcp dct iterations [%d]\n"
        "  -r, --redius <int>        span radius [%d]\n"
        "  -c, --consensus <int>     minimum consensus threshold [%d]\n"
        "  -d, --min-c-frac <float>  minimum consensus frac compared to depth [%.02f]\n"
        "  -t, --threads             thread number [%d]\n"
        "  -p, --progress            display progress\n"
        "  -v, --verbose             display messages\n"
        "  -h, --help                show this help\n"
        "\n"
        "Sketching method (default: --blend):\n"
        "      --blend               BLEND fuzzy seeds  (uses -k -w -b -n)\n"
        "      --minimizer           (w,k)-minimizers   (uses -k -w)\n"
        "      --syncmer             (s,k)-syncmers     (uses -k -s)\n"
        "      --lcp                 LCP cores          (uses -l --dct-count)\n",
        prog,
        __DEFAULT_BLEND_K__,
        __DEFAULT_BLEND_w__, 
        __DEFAULT_BLEND_B__, 
        __DEFAULT_BLEND_N__, 
        __DEFAULT_SYNCMER_S__, 
        __DEFAULT_LCP_L__, 
        __DEFAULT_LCP_D__, 
        __DEFAULT_RADIUS__, 
        __DEFAULT_CONSENSUS_THRESHOLD__, 
        __DEFAULT_CONSENSUS_FRAC_THRESHOLD__,
        __DEFAULT_THREAD_NUMBER__
    );
}

void init_params(params_t *p) {
    p->fastq = NULL;
    p->n_fastq = 0;
    p->min_consensus = __DEFAULT_CONSENSUS_THRESHOLD__;
    p->min_consensus_frac = __DEFAULT_CONSENSUS_FRAC_THRESHOLD__;
    p->radius = __DEFAULT_RADIUS__;
    p->n_threads = __DEFAULT_THREAD_NUMBER__;
    p->progress = __DEFAULT_PROGRESS__;
    p->verbose = __DEFAULT_VERBOSE__;
    // sketch
    p->sketch.type = SKETCH_BLEND;
    p->sketch.kmer_size = __DEFAULT_BLEND_K__;
    p->sketch.window_size = __DEFAULT_BLEND_w__;
    p->sketch.blend_bits = __DEFAULT_BLEND_B__;
    p->sketch.n_neighbors = __DEFAULT_BLEND_N__;
    p->sketch.smer_size = __DEFAULT_SYNCMER_S__;
    p->sketch.lcp_level = __DEFAULT_LCP_L__;
    p->sketch.dct_count = __DEFAULT_LCP_D__;
}

void parse_args(int argc, char **argv, params_t *p) {

    init_params(p);

    static struct option long_opts[] = {
        {"fasta",           required_argument, 0, 'f'},
        {"fastq",           required_argument, 0, 'q'},
        {"output",          required_argument, 0, 'o'},
        {"consensus",       required_argument, 0, 'c'},
        {"consensus-frac",  required_argument, 0, 'd'},
        {"radius",          required_argument, 0, 'r'},
        {"threads",         required_argument, 0, 't'},
        {"progress",        no_argument,       0, 'p'},
        {"verbose",         no_argument,       0, 'v'},
        // sketch
        {"kmer",            required_argument, 0, 'k'},
        {"window",          required_argument, 0, 'w'},
        {"blend-bits",      required_argument, 0, 'b'},
        {"n-neighbors",     required_argument, 0, 'n'},
        {"smer-size",       required_argument, 0, 's'},
        {"lcp-level",       required_argument, 0, 'l'},
        {"dct-count",       required_argument, 0, OPT_DCT},
        // sketch type
        {"blend",           no_argument,       0, OPT_BLEND},
        {"minimizer",       no_argument,       0, OPT_MINIMIZER},
        {"syncmer",         no_argument,       0, OPT_SYNCMER},
        {"lcp",             no_argument,       0, OPT_LCP},
        {"help",            no_argument,       0, 'h' },
        {0, 0, 0, 0}
    };

    int is_f_set = 0;
    int is_q_set = 0;
    int is_o_set = 0;
    int is_min_consensus_set = 0;
    int is_min_consensus_frac_set = 0;
    int is_type_set = 0;

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "f:q:o:k:w:b:n:s:l:r:c:d:t:pvh", long_opts, &idx)) != -1) {
        switch (opt) {
        case 'f':
            p->fasta = optarg;
            is_f_set = 1;
            break;
        case 'q': {
            /* collect all fastq files following -q: the first one (optarg)
               plus every subsequent argument until the next option or the
               end of the argument list */
            int extra = 0;
            while (optind + extra < argc && argv[optind + extra][0] != '-') {
                extra++;
            }

            /* if -q was already given, drop the previous array (last wins) */
            if (p->fastq) {
                free(p->fastq);
            }
            p->n_fastq = 1 + extra;
            p->fastq = malloc(sizeof(char *) * p->n_fastq);
            if (!p->fastq) {
                fprintf(stderr, "[%s::err] out of memory allocating fastq list\n", __TOOL_SHORT_NAME__);
                exit(1);
            }

            p->fastq[0] = optarg;
            for (int j = 0; j < extra; j++) {
                p->fastq[1 + j] = argv[optind + j];
            }
            optind += extra;   /* advance getopt past the consumed files */

            is_q_set = 1;
            break;
        }
        case 'o':
            p->bed = optarg;
            is_o_set = 1;
            break;
        case 'k':
            p->sketch.kmer_size = atoi(optarg);
            break;
        case 'w':
            p->sketch.window_size = atoi(optarg);
            break;
        case 'b':
            p->sketch.blend_bits = atoi(optarg);
            break;
        case 'n':
            p->sketch.n_neighbors = atoi(optarg);
            break;
        case 's':
            p->sketch.smer_size = atoi(optarg);
            break;
        case 'l':
            p->sketch.lcp_level = atoi(optarg);
            break;
        case OPT_DCT:
            p->sketch.dct_count = atoi(optarg);
            break;
        case 'c':
            p->min_consensus = atoi(optarg);
            is_min_consensus_set = 1;
            break;
        case 'd':
            p->min_consensus_frac = atof(optarg);
            is_min_consensus_frac_set = 1;
            break;
        case 'r':
            p->radius = atoi(optarg);
            break;
        case 't':
            p->n_threads = atoi(optarg);
            break;
        case 'p':
            p->progress = 1;
            break;
        case 'v':
            p->verbose = 1;
            break;
        case OPT_BLEND:
            p->sketch.type = SKETCH_BLEND;
            is_type_set = 1;
            break;
        case OPT_MINIMIZER:
            p->sketch.type = SKETCH_MINIMIZER;
            is_type_set = 1;
            break;
        case OPT_SYNCMER:
            p->sketch.type = SKETCH_SYNCMER;
            is_type_set = 1;
            break;
        case OPT_LCP:
            p->sketch.type = SKETCH_LCP;
            is_type_set = 1;
            break;
        case 0:
            print_usage(argv[0]);
            exit(0);
        default:
            print_usage(argv[0]);
            exit(1);
        }
    }

    (void)is_type_set;

    /* sanity checks */
    if (p->sketch.kmer_size <= 0) {
        fprintf(stderr, "[%s::err] k must be > 0 (-k <int>)\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    /* method-specific parameter validation */
    if (p->sketch.type == SKETCH_BLEND || p->sketch.type == SKETCH_MINIMIZER) {
        if (p->sketch.window_size <= 0) {
            fprintf(stderr, "[%s::err] w must be > 0 (-w <int>)\n", __TOOL_SHORT_NAME__);
            exit(1);
        }
    }

    if (p->sketch.type == SKETCH_SYNCMER) {
        if (p->sketch.smer_size <= 0) {
            fprintf(stderr, "[%s::err] s must be > 0 (-s <int>)\n", __TOOL_SHORT_NAME__);
            exit(1);
        } else if (p->sketch.smer_size >= p->sketch.kmer_size) {
            fprintf(stderr, "[%s::err] s must be < k (-s <int>)\n", __TOOL_SHORT_NAME__);
            exit(1);
        }
    }

    if (p->sketch.type == SKETCH_LCP) {
        if (p->sketch.lcp_level <= 0) {
            fprintf(stderr, "[%s::err] l must be > 0 (-l <int>)\n", __TOOL_SHORT_NAME__);
            exit(1);
        }
        if (p->sketch.dct_count < 0) {
            fprintf(stderr, "[%s::err] dct-count must be >= 0 (--dct-count <int>)\n", __TOOL_SHORT_NAME__);
            exit(1);
        }
    }

    if (!is_f_set) {
        fprintf(stderr, "[%s::err] fasta is not provided (-f <file>)\n", __TOOL_SHORT_NAME__);
        print_usage(argv[0]);
        exit(1);
    }

    if (!is_q_set) {
        fprintf(stderr, "[%s::err] query is not provided (-q <file>)\n", __TOOL_SHORT_NAME__);
        print_usage(argv[0]);
        exit(1);
    }

    if (!is_o_set) {
        fprintf(stderr, "[%s::err] output file is not provided (-o <file>)\n", __TOOL_SHORT_NAME__);
        print_usage(argv[0]);
        exit(1);
    }

    if (!file_exists(p->fasta)) {
        fprintf(stderr, "[%s::err] fasta file does not exists (-f <file>)\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    for (int i = 0; i < p->n_fastq; i++) {
        if (!file_exists(p->fastq[i])) {
            fprintf(stderr, "[%s::err] fastq file does not exists (%s)\n", __TOOL_SHORT_NAME__, p->fastq[i]);
            exit(1);
        }
    }

    if (p->n_threads <= 0 || 1024 < p->n_threads) {
        fprintf(stderr, "[%s::err] invalid thread count (-t <int>)\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    if (is_min_consensus_frac_set) {
        p->use_consensus_frac = 1;
    } else if (is_min_consensus_set) {
        p->use_consensus_frac = 0;
    } else {
        p->use_consensus_frac = 1;
    }

    // derive runtime layout (kmer_shift / span_mask) from the selected method
    sketch_args_finalize(&p->sketch);

    if (p->verbose) {
        fprintf(stderr, "[%s::args] fa: %s\n", __TOOL_SHORT_NAME__, p->fasta);
        for (int i = 0; i < p->n_fastq; i++) {
            fprintf(stderr, "[%s::args] fq[%d]: %s\n", __TOOL_SHORT_NAME__, i, p->fastq[i]);
        }
        fprintf(stderr, "[%s::args] method: %s\n", __TOOL_SHORT_NAME__, sketch_type_name(p->sketch.type));
        switch (p->sketch.type) {
        case SKETCH_BLEND:
            fprintf(stderr, "[%s::args] k: %d, w: %d, b: %d, n: %d\n",
                    __TOOL_SHORT_NAME__, p->sketch.kmer_size, p->sketch.window_size, p->sketch.blend_bits, p->sketch.n_neighbors);
            break;
        case SKETCH_MINIMIZER:
            fprintf(stderr, "[%s::args] k: %d, w: %d\n",
                    __TOOL_SHORT_NAME__, p->sketch.kmer_size, p->sketch.window_size);
            break;
        case SKETCH_SYNCMER:
            fprintf(stderr, "[%s::args] k: %d, s: %d\n",
                    __TOOL_SHORT_NAME__, p->sketch.kmer_size, p->sketch.smer_size);
            break;
        case SKETCH_LCP:
            fprintf(stderr, "[%s::args] l: %d, dct: %d\n",
                    __TOOL_SHORT_NAME__, p->sketch.lcp_level, p->sketch.dct_count);
            break;
        default:
            break;
        }
        fprintf(stderr,
            "[%s::args] r: %d, %c: %0.2f, t: %d\n",
            __TOOL_SHORT_NAME__, p->radius,
            p->use_consensus_frac ? 'd' : 'c',
            p->use_consensus_frac ? p->min_consensus_frac : p->min_consensus,
            p->n_threads
        );
    }
}
