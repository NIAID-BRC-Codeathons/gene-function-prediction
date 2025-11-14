/*  Read codon counts of a genome.  Then read pairs of codon usage frequencies
 *  from stdin, reporting the axis projection position and P-value for each
 *  gene.
 *
 *     project_codon_usage_onto_axis [-l max_len] codon_counts < codon_freq_pairs > x_tab_p
 *
 *  Options:
 *
 *    -l max_len   #  scale chisquare *= max_len / n_codons for large genes
 *    -v           #  print version number and exit
 *
 *  Count file can be produced with compile_codon_counts:
 *
 *     compile_codon_counts < cds_fasta_file > codon_count_file
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gjocodonlib.c"  /* This will bring in gjocodonlib.h */

/* Types */

typedef struct { double chisqr; int df; } chisqr_df;

/* Prototypes */

int   read_cds_codon_counts( FILE *fp, codon_cnt_t cds[], int nmax );
void  project_on_freq_pairs( FILE *fp, int ncds, codon_cnt_t cds[], double max_l );
void  usage( char * command );

/* Global values */

#define  VERSION  "1.00"
#define  MAX_CDS   64000
#define  MIN_FREQ  0.00001

codon_cnt_t  cds[MAX_CDS];  /*  Array of codon counts  */

int  n_aa = 18;


int
main( int argc, char** argv )
{
    FILE *fp;
    int    i, ncds;
    double max_l =  0.0;  /*  No length limit in scoring  */

    /* Read command flags */

    i = 1;
    while ( ( i < argc ) && ( argv[i][0] == '-' ) )
    {
        if ( argv[i][1] == 'l' )
        {
            if ( argv[i][2] ) { argv[i] += 2; }
            else if ( ++i >= argc )
            {
                fprintf( stderr, "Missing value after -l\n" );
                usage( argv[0] );
                return 0;
            }
            max_l = atof( argv[i] );
            // fprintf( stderr, "Maximum length set to %.0f.\n", max_l );
        }
        else if ( argv[i][1] == 'v' )
        {
            printf( "%s\n", VERSION );
            return 1;
        }
        else
        {
            fprintf( stderr, "Bad flag '%s'.\n", argv[i] );
            usage( argv[0] );
            return 0;
        }

        i++;
    }

    if ( argc < i+1 )
    {
        fprintf( stderr, "Missing codon count file name.\n" );
        usage( argv[0] );
        return 0;
    }

    fp = fopen( argv[i], "r" );
    if ( ! fp )
    {
        fprintf( stderr, "Could not open codon count file '%s'.\n", argv[i] );
        usage( argv[0] );
        return 0;
    }

    if ( argc > i+1 )
    {
        fprintf( stderr, "Extra parameter after codon count file: '%s'.\n", argv[i+1] );
        usage( argv[0] );
        return 0;
    }

    ncds = read_cds_codon_counts( fp, cds, MAX_CDS );
    ( void ) fclose( fp );

    // fprintf( stderr, "Number of codon counts (genes) read = %d\n", ncds );

    if ( ! ncds )
    {
        fprintf( stderr, "Failed to read codon counts from '%s'.\n", argv[i] );
        usage( argv[0] );
        return 0;
    }

    project_on_freq_pairs( stdin, ncds, cds, max_l );
    return 1;
}


/*-----------------------------------------------------------------------------
 *  Read the file of white-space separated codon counts for the genome:
 *
 *  int read_cds_codon_counts( FILE *fp, codon_cnt_t cds[], int nmax )
 *
 *-----------------------------------------------------------------------------
 */

int
read_cds_codon_counts( FILE *fp, codon_cnt_t cds[], int nmax )
{
    int  ncds = 0;

    while ( ncds < nmax )
    {
        cds[ncds] = read_codon_counts( fp );
        if ( cds[ncds].naa < n_aa ) { break; }
        ncds++;
    }

    return ncds;
}


void
project_on_freq_pairs( FILE *fp, int ncds, codon_cnt_t cds[], double max_l )
{
    codon_freq_t  freq1, freq2;
    k_vector_t    k;
    x_pval_t      x_p;
    int           i;
    int           ok = 1;

    while ( ok )
    {
        freq1 = read_codon_freqs( fp );
        if ( freq1.naa < 18 ) { break; }
        freq2 = read_codon_freqs( fp );
        if ( freq2.naa < 18 ) { break; }

        k = k_from_f0_and_f1( freq1, freq2 );

        for ( i = 0; i < ncds; i++ )
        {
            x_p = project_by_min_chisqr( freq1, k, cds[i], max_l );
            printf( "%.3f\t%.5e\n", x_p.x, x_p.p );
        }
        fflush( stdout );
    }
}


/*-----------------------------------------------------------------------------
 *  Usage statement:
 *
 *      void usage( char * command )
 *
 *-----------------------------------------------------------------------------
 */

void
usage( char * command )
{
    fprintf( stderr, "\n" );
    fprintf( stderr, "Usage: %s [options] codon_count_file < codon_frequency_pairs > x_tab_p\n", command );
    fprintf( stderr, "\n" );
    fprintf( stderr, "Options:\n\n" );
    fprintf( stderr, "   -l max_len   #  In calculating P-value, limit codons to max_len\n" );
    fprintf( stderr, "   -v           #  Print version number and exit\n" );
    fprintf( stderr, "\n" );

    return;
}

