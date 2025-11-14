/*  Read codon counts of a genome.  Then read pairs of codon usage frequencies
 *  from stdin, reporting the axis projection position and P-value for each
 *  gene.
 *
 *     codon_counts_x_and_p [-l max_len] < freq_pairs_and_codon_counts > x_tab_p
 *
 *  Options:
 *
 *    -d           #  include id and definition
 *    -i           #  include id
 *    -l max_len   #  scale chisquare *= max_len / n_codons for large genes
 *    -v           #  print version number and exit
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "gjocodonlib.c"  /* This will bring in gjocodonlib.h */

/* Types */

typedef struct { double chisqr; int df; } chisqr_df;

/* Prototypes */

int   read_cds_codon_counts( FILE *fp, codon_cnt_t cds[], int nmax );
void  project_on_freq_pairs( FILE *fp, int ncds, codon_cnt_t cds[], double max_l );
void  usage( char * command );

/* Global values */

#define  VERSION  "1.01"
#define  MIN_FREQ  0.00001
#define  MAX_LABEL_LEN  4999

int  n_aa = 18;


int
main( int argc, char** argv )
{
    codon_freq_t  freq1, freq2;
    k_vector_t    k;
    codon_cnt_t   cds;  /*  Codon counts  */
    x_pval_t      x_p;
    double        max_l = 0.0;  /*  No length limit in scoring  */
    int           ok    = 1;
    int           id    = 0;
    int           i;
    char          label[MAX_LABEL_LEN+1];


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
        else if ( argv[i][1] == 'd' )   // include the definition
        {
            id = 2;
        }
        else if ( argv[i][1] == 'i' )   // include the id
        {
            id = 1;
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

    /* Read 2 frequencies */

    freq1 = read_codon_freqs( stdin );
    if ( freq1.naa < n_aa ) { ok = 0; }
    else
    {
        freq2 = read_codon_freqs( stdin );
        if ( freq2.naa < n_aa ) { ok = 0; }
    }

    if ( ! ok )
    {
        fprintf( stderr, "Failed to read codon frequencies.\n" );
        usage( argv[0] );
        return 0;
    }

    k = k_from_f0_and_f1( freq1, freq2 );

    while ( ok )
    {
        cds = read_codon_counts_2( stdin, label, MAX_LABEL_LEN );
        if ( cds.naa < n_aa ) { break; }
        x_p = project_by_min_chisqr( freq1, k, cds, max_l );
        if ( id && label[0] )
        {
            if ( id == 1 )  //  For id only, terminate at white space
            {
                char * cp = label;
                while ( *cp && ! isspace( *cp ) ) { cp++; }
                *cp++ = '\n';
                *cp   = '\0';
            }
            printf( "%.3f\t%.5e\t%s", x_p.x, x_p.p, label );
        }
        else
        {
            printf( "%.3f\t%.5e\n", x_p.x, x_p.p );
        }
        fflush( stdout );
    }

    return 1;
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
    fprintf( stderr, "Usage: %s [options] < frequency_pairs_cdon_counts > x_tab_p\n", command );
    fprintf( stderr, "\n" );
    fprintf( stderr, "Options:\n\n" );
    fprintf( stderr, "   -d           #  Include the definition\n" );
    fprintf( stderr, "   -i           #  Include the id\n" );
    fprintf( stderr, "   -l max_len   #  In calculating P-value, limit codons to max_len\n" );
    fprintf( stderr, "   -v           #  Print version number and exit\n" );
    fprintf( stderr, "\n" );

    return;
}

