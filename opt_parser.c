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
        "  -k, --kmer <int>          k-mer size [%d]\n"
        "  -w, --window <int>        window size [%d]\n"
        "  -b, --blend-bits <int>    number of bits for hash [%d]\n"
        "  -n, --n-neighbors <int>   number of neighbour [%d]\n"
        "  -c, --consensus <int>     minimum consensus threshold [%d]\n"
        "  -t, --threads             thread number [%d]\n"
        "  -p, --progress            display progress\n"
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
    p->n_threads = __DEFAULT_THREAD_NUMBER__;
}

void parse_args(int argc, char **argv, params_t *p) {

    init_params(p);

    static struct option long_opts[] = {
        {"fasta",       required_argument, 0, 'f'},
        {"fastq",       required_argument, 0, 'q'},
        {"kmer",        required_argument, 0, 'k'},
        {"window",      required_argument, 0, 'w'},
        {"blend-bits",  required_argument, 0, 'b'},
        {"n-neighbors", required_argument, 0, 'n'},
        {"consensus",   required_argument, 0, 'c'},
        {"threads",     required_argument, 0, 't'},
        {"progress",    no_argument,       0, 'p'},
        {"help",        no_argument,       0, 'h' },
        {0, 0, 0, 0}
    };

    int is_f_set = 0; 
    int is_q_set = 0;

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "f:q:k:w:b:n:c:t:ph", long_opts, &idx)) != -1) {
        switch (opt) {
        case 'f':
            p->fasta = optarg;
            is_f_set = 1;
            break;
        case 'q':
            p->fastq = optarg;
            is_q_set = 1;
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
        fprintf(stderr, "[ERROR] k and w must be > 0 (-k <int>)\n");
        exit(1);
    }

    if (p->w <= 0) {
        fprintf(stderr, "[ERROR] w must be > 0 (-w <int>)\n");
        exit(1);
    }

    if (!is_f_set) {
        fprintf(stderr, "[ERROR] fasta is not provided (-f <file>)\n");
        print_usage(argv[0]);
        exit(1);
    }

    if (!is_q_set) {
        fprintf(stderr, "[ERROR] query is not provided (-q <file>)\n");
        print_usage(argv[0]);
        exit(1);
    }

    if (!file_exists(p->fasta)) {
        fprintf(stderr, "[ERROR] fasta file does not exists (-f <file>)\n");
        exit(1);
    }

    if (!file_exists(p->fastq)) {
        fprintf(stderr, "[ERROR] fastq file does not exists (-q <file>)\n");
        exit(1);
    }

    if (p->n_threads <= 0 || 1024 < p->n_threads) {
        fprintf(stderr, "[ERROR] invalid thread count (-t <int>)\n");
        exit(1);  
    }
}