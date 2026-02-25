#include "opt_parser.h"


int file_exists(char *filename) {
    struct stat buffer;   
    return (stat(filename, &buffer) == 0);
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
        "  -c, --consensus <int>     minimum consensus threshold [%d]\n"
        "  -t, --threads             thread number [%d]\n"
        "  -p, --progress            display progress\n"
        "  -v, --verbose             display messages\n"
        "  -h, --help                show this help\n",
        prog, __DEFAULT_BLEND_K__, __DEFAULT_BLEND_w__, __DEFAULT_BLEND_BITS__, __DEFAULT_BLEND_NEIGHBOR_NUMBER__, __DEFAULT_CONSENSUS_THRESHOLD__, __DEFAULT_THREAD_NUMBER__
    );
}

void init_params(params_t *p) {
    p->k = __DEFAULT_BLEND_K__;
    p->w = __DEFAULT_BLEND_w__;
    p->blend_bits = __DEFAULT_BLEND_BITS__;
    p->n_neighbors = __DEFAULT_BLEND_NEIGHBOR_NUMBER__;
    p->min_consensus = __DEFAULT_CONSENSUS_THRESHOLD__;
    p->progress = __DEFAULT_PROGRESS__;
    p->verbose = __DEFAULT_VERBOSE__;
    p->n_threads = __DEFAULT_THREAD_NUMBER__;
}

void parse_args(int argc, char **argv, params_t *p) {

    init_params(p);

    static struct option long_opts[] = {
        {"fasta",       required_argument, 0, 'f'},
        {"fastq",       required_argument, 0, 'q'},
        {"output",      required_argument, 0, 'o'},
        {"kmer",        required_argument, 0, 'k'},
        {"window",      required_argument, 0, 'w'},
        {"blend-bits",  required_argument, 0, 'b'},
        {"n-neighbors", required_argument, 0, 'n'},
        {"consensus",   required_argument, 0, 'c'},
        {"threads",     required_argument, 0, 't'},
        {"progress",    no_argument,       0, 'p'},
        {"verbose",     no_argument,       0, 'v'},
        {"help",        no_argument,       0, 'h' },
        {0, 0, 0, 0}
    };

    int is_f_set = 0; 
    int is_q_set = 0;
    int is_o_set = 0;

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "f:q:o:k:w:b:n:c:t:pvh", long_opts, &idx)) != -1) {
        switch (opt) {
        case 'f':
            p->fasta = optarg;
            is_f_set = 1;
            break;
        case 'q':
            p->fastq = optarg;
            is_q_set = 1;
            break;
        case 'o':
            p->bed = optarg;
            is_o_set = 1;
            break;
        case 'k':
            p->k = atoi(optarg);
            break;
        case 'w':
            p->w = atoi(optarg);
            break;
        case 'b':
            p->blend_bits = atoi(optarg);
            break;
        case 'n':
            p->n_neighbors = atoi(optarg);
            break;
        case 'c':
            p->min_consensus = atoi(optarg);
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
        case 0:
            print_usage(argv[0]);
            exit(0);
        default:
            print_usage(argv[0]);
            exit(1);
        }
    }

    /* sanity checks */
    if (p->k <= 0) {
        fprintf(stderr, "[%s::err] k and w must be > 0 (-k <int>)\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    if (p->w <= 0) {
        fprintf(stderr, "[%s::err] w must be > 0 (-w <int>)\n", __TOOL_SHORT_NAME__);
        exit(1);
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

    if (!file_exists(p->fastq)) {
        fprintf(stderr, "[%s::err] fastq file does not exists (-q <file>)\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    if (p->n_threads <= 0 || 1024 < p->n_threads) {
        fprintf(stderr, "[%s::err] invalid thread count (-t <int>)\n", __TOOL_SHORT_NAME__);
        exit(1);  
    }

    if (p->verbose) {
        printf(
            "[%s:args] fa: %s\n"
            "[%s:args] fq: %s\n"
            "[%s:args] k: %d, w: %d, b: %d, n: %d, c: %d, t: %d\n",
            __TOOL_SHORT_NAME__, p->fasta, 
            __TOOL_SHORT_NAME__, p->fastq, 
            __TOOL_SHORT_NAME__, p->k, p->w, p->blend_bits, p->n_neighbors, p->min_consensus, p->n_threads
        );
    }
}