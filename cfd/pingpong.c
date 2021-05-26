#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>

#define PING 1
#define PONG 2

static void print_usage(void);
static void print_version(void);
static void print_help(void);

static char *progname;

#define PACKAGE "pingpong"
#define VERSION "1.0"

/* Command line options */
static struct option long_opts[] = {
    { "help",    0, NULL, 'h' },
    { "version", 0, NULL, 'V' },
    { "minsize", 1, NULL, 'm' },
    { "maxsize", 1, NULL, 'n' },
    { "count",   1, NULL, 'c' },
    { 0,         0, 0,    0   } 
};

#define GETOPTS "c:hm:n:V"

int main(int argc, char **argv)
{
    int i, n, p, size;
    int iters = 1000, minsize = 1, maxsize = 32768;
    int show_help = 0, show_usage = 0, show_version = 0;
    double start;

    progname = argv[0];
    
    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch (optc) {
            case 'h':
                show_help = 1;
                break;
            case 'V':
                show_version = 1;
                break;
            case 'm':
                minsize = atoi(optarg);
                if (minsize < 1 || minsize > 65536) {
                    show_usage = 1;
                }
                break;
            case 'n':
                maxsize = atoi(optarg);
                if (maxsize < 1 || maxsize > 65536) {
                    show_usage = 1;
                }
                break;
            case 'c':
                printf("iters: %d\n", iters);
                if (iters < 1) {
                    show_usage = 1;
                }
                break;
            default:
                show_usage = 1;
        }
    }

    if (show_version) {
        print_version();
        if (!show_help) {
            return 0;
        }
    }

    if (show_help) {
        print_help();
        return 0;
    }

    if (show_usage || optind < argc || maxsize < minsize) {
        print_usage();
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);


    if (n != 2) {
        fprintf(stderr, "%s must be run on exactly 2 processors.\n", progname);
        MPI_Finalize();
        return 1;
    }

    char *dummy;
    dummy = malloc(maxsize);
    for (size = minsize; size <= maxsize; size *= 2) {
        start = MPI_Wtime();
        for (i = 0; i < iters; i++) {
            MPI_Status s;
            /* ping... */
            if (p == 0) {
                MPI_Send(dummy, size, MPI_CHAR, 1, PING, MPI_COMM_WORLD);
            } else {
                MPI_Recv(dummy, size, MPI_CHAR, 0, PING, MPI_COMM_WORLD, &s);
            }

            /* pong... */
            if (p == 1) {
                MPI_Send(dummy, size, MPI_CHAR, 0, PONG, MPI_COMM_WORLD);
            } else {
                MPI_Recv(dummy, size, MPI_CHAR, 1, PONG, MPI_COMM_WORLD, &s);
            }
        }
        if (p == 0) { 
            printf("%d %d byte pingpongs took %f seconds.\n", iters,
                size, MPI_Wtime() - start);
        }
    }
    MPI_Finalize();
    return 0;
}

static void print_usage(void)
{
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
}

static void print_version(void)
{
    fprintf(stderr, "%s %s\n", PACKAGE, VERSION);
}

static void print_help(void)
{
    fprintf(stderr, "%s. A utility for benchmarking MPI communications.\n\n",
        PACKAGE);
    fprintf(stderr, "Usage %s [OPTIONS]\n\n", progname);
    fprintf(stderr, "  -h, --help           Print a summary of the options\n");
    fprintf(stderr, "  -V, --version        Print the version number\n");
    fprintf(stderr, "  -m, --minsize        The minimum packet size. Must be 1-65536\n");
    fprintf(stderr, "  -n, --maxsize        The maximum packet size. Must be 1-65536\n");
    fprintf(stderr, "  -c, --count          The number of pings to send for a given size\n");
}
