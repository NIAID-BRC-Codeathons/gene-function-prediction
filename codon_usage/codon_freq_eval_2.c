/*  Read codon counts of a genome.  Then evaluate the agreement of the
 *  genome genes with codon usage frequencies read from stdin.  If a P-value
 *  is not supplied, each score is the sum of the chi-square P-values for the
 *  set of genes matching the frequencies.  If a P-value is supplied, each
 *  score is the number of genes matching the frequencies with P >= P-value.
 *
 *     codon_freq_eval_2 [-e exp] [-l max_len] [-p pval] codon_count_file < codon_frequencies > scores
 *
 *  Options:
 *
 *    -e expon     #  score for each codon count is p_value**expon
 *    -l max_len   #  scale chisquare *= max_len / n_codons for large genes
 *    -p p_value   #  score is total genes with P-value >= p_value
 *    -v           #  print version number and exit
 *
 *  Count file can be produced with compile_codon_counts:
 *
 *     compile_codon_counts < cds_fasta_file > codon_count_file
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Types */

typedef struct { double chisqr; int df; } chisqr_df;

/* Prototypes */

int        read_cds_codon_counts( FILE *fp, int cds[][64], int len[] );
void       evaluate_codon_frequencies( int ncds, int cds[][64], int *len, double p_val, double expon, double max_l );
double     evaluate_freq( double *freq, int ncds, int cds[][64], int *len, double p_val, double expon, double max_l );
void       usage( char * command );
chisqr_df  chi_square( int n, double * expected, int * counts );
double     chisqr_prob( double chisqr, int df );

/* Global values */

#define  version  "1.00"
#define  max_cds   64000
#define  min_freq  0.000001

int  cds[max_cds][64];  /*  Array of codon counts  */
int  len[max_cds];      /*  Array of lengths  */

                   /* A  C  D  E  F  G  H  I  K  L  N  P  Q  R  S  T  V  Y */
int  aa_codon[18] = { 4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 2, 4, 2, 6, 6, 4, 4, 2 };
int  n_aa    = 18;
int  n_codon = 59;


int main( int argc, char** argv )
{
    FILE *fp;
    int    i, ncds;
    double p_val = -1.0;  /*  Do not score by p-value threshold  */
    double expon =  1.0;  /*  Sum of p ** expon  */
    double max_l =  0.0;  /*  No length limit in scoring  */

    /* Read command flags */

    i = 1;
    while ( ( i < argc ) && ( argv[i][0] == '-' ) )
    {
        if ( argv[i][1] == 'e' )
        {
            if ( argv[i][2] ) { argv[i] += 2; }
            else if ( ++i >= argc )
            {
                fprintf( stderr, "Missing value after -e\n" );
                usage( argv[0] );
                return 0;
            }
            expon = atof( argv[i] );
            // fprintf( stderr, "Exponent set to %.2e.\n", expon );
        }
        else if ( argv[i][1] == 'l' )
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
        else if ( argv[i][1] == 'p' )
        {
            if ( argv[i][2] ) { argv[i] += 2; }
            else if ( ++i >= argc )
            {
                fprintf( stderr, "Missing value after -p\n" );
                usage( argv[0] );
                return 0;
            }
            p_val = atof( argv[i] );
            // fprintf( stderr, "P-value set to %.1e.\n", p_val );
        }
        else if ( argv[i][1] == 'v' )
        {
            printf( "%s\n", version );
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

    ncds = read_cds_codon_counts( fp, cds, len );
    ( void ) fclose( fp );

    // fprintf( stderr, "Number of codon counts (genes) read = %d\n", ncds ); return 1;

    if ( ! ncds )
    {
        fprintf( stderr, "Failed to read codon counts from '%s'.\n", argv[i] );
        usage( argv[0] );
        return 0;
    }

    evaluate_codon_frequencies( ncds, cds, len, p_val, expon, max_l );
    return 1;
}


/*-----------------------------------------------------------------------------
 *  Read the file of white-space separated codon counts for the genome:
 *
 *      int  read_cds_codon_counts( FILE *fp, int cds[][64], int *len )
 *
 *-----------------------------------------------------------------------------
 */

int  read_cds_codon_counts( FILE *fp, int cds[][64], int *len )
{
    int  codons_read, c;
    int *ptr;
    int  ok   = 1;
    int  ncds = 0;
    int  n;

    while ( ok == 1 )
    {
        ptr = &(cds[ncds][0]);
        codons_read = 0;
        n = 0;
        while ( ( ok == 1 ) && ( codons_read < n_codon ) )
        {
            ok = fscanf( fp, "%i", ptr );
            if ( ok == 1 )
            {
                n += *ptr;
                codons_read++;
                ptr++;
            }
        }

        if ( ok == 1 )
        {
            *len++ = n;   /*  Record the protein length  */

            /*  Flush to end of line so that extra codons or an id are okay  */

            while ( ( ( c = getc( fp ) ) > 0 ) && ( c != '\n' ) ) {}
            if ( ( ++ncds >= max_cds ) || ( c < 0 ) ) { ok = 0; }
        }
    }

    return ncds;
}


void  evaluate_codon_frequencies( int ncds, int cds[][64], int *len,
                                  double p_val, double expon, double max_l
                                )
{
    double  freq[64];
    double *inptr;
    double  score;
    int  codons_read;
    int  i, c;
    int  ok = 1;

    while ( ok )
    {
        inptr = &(freq[0]);
        codons_read = 0;

        /* check for EOF or blank line */
        c = getc( stdin );
        if ( c < 0 ) { break; }

        ungetc( c, stdin );
        if ( c == '\n' )
        {
            score = -1;
        }
        else
        {
            while ( ok && ( codons_read < n_codon ) )
            {
                ok = scanf( "%lf", inptr );
                if ( ! ok ) { break; }
                if ( *inptr < min_freq ) { *inptr = min_freq; }
                codons_read++;
                inptr++;
            }
            if ( ok )
            {
                score = evaluate_freq( freq, ncds, cds, len, p_val, expon, max_l );
            }
            else
            {
                score = -1;
            }
        }

        fprintf( stdout, "%.8f\n", score );
        fflush( stdout );

        if ( ! ok ) { break; }

        /*  Flush the input line  */
        while ( ( ( c = getc( stdin ) ) > 0 ) && ( c != '\n' ) ) {}
        if ( c < 0 ) { break; }
    }
}


double  evaluate_freq( double *freq, int ncds, int cds[][64], int *len,
                       double p_val, double expon, double max_l
                     )
{
    chisqr_df  chi2;
    double     p, ttl_scr, ttl_chisqr, *freqptr;
    int       *cdsptr;
    int        aa, n, i, df;

    ttl_scr = 0.0;
    for ( i = 0; i < ncds; i++ )
    {
        freqptr = freq;
        cdsptr  = &(cds[i][0]);
        df = 0;
        ttl_chisqr = 0.0;
        for ( aa = 0; aa < n_aa; aa++ )
        {
            n = aa_codon[ aa ];
            chi2 = chi_square( n, freqptr, cdsptr );
            if ( chi2.df > 0 && chi2.chisqr >= 0 )
            {
                ttl_chisqr += chi2.chisqr;
                df         += chi2.df;
            }
            freqptr += n;
            cdsptr  += n;
        }
        if ( df > 0 )
        {
            /*  Adjust effective length, if necessary  */
            if ( ( max_l > 0 ) && ( len[i] > max_l ) )
            {
                ttl_chisqr *= ( max_l / len[i] );
            }

            /*  Evaluate probability of chi-square  */
            p = chisqr_prob( ttl_chisqr, df );

            /*  Record desired score  */
            ttl_scr += p_val < 0  ? ( expon == 1 ? p : pow( p, expon ) )
                     : p >= p_val ? 1
                     :              0;
        }
    }

    return ttl_scr;
}


/*-----------------------------------------------------------------------------
 *  Usage statement:
 *
 *      void usage( char * command )
 *
 *-----------------------------------------------------------------------------
 */

void usage( char * command )
{
    fprintf( stderr, "\n" );
    fprintf( stderr, "Usage: %s [options] codon_count_file < codon_frequencies > scores\n", command );
    fprintf( stderr, "\n" );
    fprintf( stderr, "Options:\n\n" );
    fprintf( stderr, "   -e exponent  #  Score is sum over genes of p ** exponent (D = 1)\n" );
    fprintf( stderr, "   -l max_len   #  In calculating P-value, limit gene length to max_len\n" );
    fprintf( stderr, "   -p p_value   #  Score is the count of genes with p >= p_value\n" );
    fprintf( stderr, "\n" );

    return;
}


/*-----------------------------------------------------------------------------
 *  Find the chi-square value for an expected value (or frequency) list and
 *  an observed list (expected values need not be normalized):
 *
 *      double chi_sqr = chi_square( n, expected, counts )
 *
 *  Negative degrees of freedom signifies error.
 *-----------------------------------------------------------------------------
 */

chisqr_df  chi_square( int n, double * expected, int * counts )
{
    chisqr_df  result = { 0.0,  0 };
    chisqr_df  error  = { 0.0, -1 };
    double  total_exp = 0;
    double  e, e_scale, chisqr, diff;
    int     df;
    int     total_cnt = 0;
    int     c, i;

    if ( n < 2 )                  { return result; }
    if ( ! expected || ! counts ) { return error; }

    df = -1;
    for ( i = 0; i < n; i++ )
    {
        e = expected[i];
        if ( e < 0 ) { return error; }
        c = counts[i];
        if ( ( c < 0 ) || ( e == 0 && c > 0 ) ) { return error; }

        if ( e > 0 )
        {
            total_exp += e;
            total_cnt += c;
            df++;
        }
    }

    if ( df <= 0 || total_cnt == 0 ) { return result; }  /* no df */
    if ( total_exp == 0 )            { return error; }

    e_scale = total_cnt / total_exp;
    chisqr = 0;
    for ( i = 0; i < n; i++ )
    {
        e = *expected++ * e_scale;
        if ( e > 0 )
        {
            diff = *counts - e;
            chisqr += diff * diff / e;
        }
        counts++;      /* need to consume counts regardless of e */
    }

    result.chisqr = chisqr;
    result.df     = df;
    return result;
}


/*-----------------------------------------------------------------------------
 *   Probability of a chi square value greater than or equal to chisqr
 *
 *   Based on:
 *
 *   Zelen, M. and Severo, N. C. (1965).  Probability functions.  In
 *   Handbook of Mathematical Functions, Abramowitz, M. and Stegun,
 *   I. A., eds. (New York: Dover Publications), pp. 925-995.
 *
 *   Programmed in C by Gary Olsen
 *
 *      double p_value = chisqr_prob( double chisqr, int df )
 *
 *-----------------------------------------------------------------------------
 */

double chisqr_prob( double chisqr, int df )
{
    double  sum, delta, denom1, denom2, prob_not;
    int     i;

    if ( ( chisqr < 0 ) || ( df <= 0 ) ) { return -1; }
    if ( chisqr == 0 ) { return 1; }
    if ( chisqr - df - 10*sqrt(df) > 49.0 ) { return 1e-14; }

    sum    = 1;
    delta  = 1;
    denom1 = df;

    while ( 1e16 * delta > sum ) {
        denom1 += 2;
        delta  *= chisqr / denom1;
        sum    += delta;
    }

    denom2 = df;
    for ( i = df - 2; i > 0; i -= 2 ) { denom2 *= i; }

    prob_not = sum * exp(-0.5 * chisqr) * pow( chisqr, (df+1)/2 ) / denom2;

    /*  pi/2 = 1.57079632679489661923  */

    if ( ( df % 2 ) != 0 ) { prob_not /= sqrt( 1.57079632679489661923 * chisqr ); }

    return 1.0 - prob_not;
}

