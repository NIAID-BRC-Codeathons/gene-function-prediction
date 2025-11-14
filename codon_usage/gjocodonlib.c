#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "gjocodonlib.h"

/*==============================================================================
 *  Summary of functions (prototypes in gjocodonlib.h)
 *
 *  double        count_vs_freq_pval( codon_cnt_t cnt, codon_freq_t freq, int maxlen );
 *  chisqr_n_t    count_vs_freq_chisqr( codon_cnt_t cnt, codon_freq_t freq );
 *  
 *  x_pval_t      project_by_min_chisqr( codon_freq_t f0, k_vector_t k, codon_cnt_t counts, int maxlen );
 *  codon_freq_t  freqs_at_x( codon_freq_t f0, k_vector_t k, double x );
 *  k_vector_t    k_from_f0_and_f1( codon_freq_t f0, codon_freq_t f1 );
 *  double        pval_at_x( codon_cnt_t counts, int maxlen, codon_freq_t f0, k_vector_t k, double x );
 *  chisqr_n_t    chisqr_at_x( codon_cnt_t counts, codon_freq_t f0, k_vector_t k, double x );
 *  codon_freq_t  freqs_at_x_given_f0_and_f1( codon_freq_t f0, codon_freq_t f1, double x );
 *  
 *  codon_cnt_t   read_codon_counts( FILE * fp );
 *  codon_cnt_t   parse_codon_counts( char * cp );
 *  int           s_next_count( char ** cpp );
 *  int           s_get_digit( char ** cpp );
 *  int           s_flush_count_space( char ** cpp );
 *  
 *  codon_freq_t  read_codon_freqs( FILE * fp );
 *  codon_freq_t  parse_codon_freqs( char * cp );
 *  double        s_next_frequency( char ** cpp );
 *  int           s_flush_freq_div( char ** cpp );
 *  
 *  void          fprint_codon_counts( FILE * fp, codon_cnt_t counts );
 *  void          fprint_codon_freqs( FILE * fp, codon_freq_t freqs );
 *  
 *  int           set_codons_per_aa( int naa, int * new_aa_n_cdn );
 *  int           reset_codons_per_aa();
 *  
 *  chisqr_n_t    chisqr( double * freq, int * cnt, int n );
 *  double        chisqr_prob( double chisqr, int df );
 *  double        chisqr_prob_w_maxn( double chi2, int df, int n, int maxn );
 *  
 *==============================================================================
 *  Basic p-value calculation
 *-----------------------------------------------------------------------------
 *  Compare codon counts to expected frequencies by chi-square.
 *
 *     pval = count_vs_freq_pval( cncodon_cnt_t cnt, codon_freq_t freq, int maxlen )
 *
 *-----------------------------------------------------------------------------
 */
double
count_vs_freq_pval( codon_cnt_t cnt, codon_freq_t freq, int maxlen )
{
    chisqr_n_t  chi2_n;
    chi2_n = count_vs_freq_chisqr( cnt, freq );
    return chisqr_prob_w_maxn( chi2_n.chisqr, chi2_n.df, chi2_n.n, maxlen );
}


/*-----------------------------------------------------------------------------
 *  Compare codon counts to expected frequencies by chi-square.
 *
 *     ( $chisqr, $df, $n ) = count_vs_freq_chisqr( codon_cnt_t cnt, codon_freq_t freq )
 *
 *-----------------------------------------------------------------------------
 */
chisqr_n_t
count_vs_freq_chisqr( codon_cnt_t cnt, codon_freq_t freq )
{
    chisqr_n_t  result;
    chisqr_n_t  partial;
    double      sum = 0;
    int         df = 0;
    int         total = 0;
    int         i, offset;

    for ( i = 0; i <= 17; i++ )
    {
        offset = aa_offset[i];
        partial = chisqr( &freq.freq[offset], &cnt.cnt[offset], aa_n_cdn[i] );
        if ( partial.df > 0 )
        {
            sum   += partial.chisqr;
            df    += partial.df;
            total += partial.n;
        }
    }

    result.chisqr = sum;
    result.df     = df;
    result.n      = total;

    return result;
}


/*==============================================================================
 *  Projections of counts onto an axis defined by two frequencies.
 *-----------------------------------------------------------------------------
 *  Description from the perl version:
 *
 *  For each frequency in the arguements, find the position along the
 *  $freq0 to $freq1 vector that gives the smallest chi square value.
 *  Coordinates along vector are measured with $freq0 being 0, and $freq1
 *  being 1.
 *
 *     @projections = project_on_freq_vector_by_chi_sqr_2( \@f0, \@f1,   \@cnts1, \@cnts2, ...   )
 *     @projections = project_on_freq_vector_by_chi_sqr_2( \@f0, \@f1, [ \@cnts1, \@cnts2, ... ] )
 *
 *  The returned projections are quadruples or pentuples:
 *
 *     [ $projection_on_f0_f1_axis, $chi_square, $deg_of_freedom, $n_codon ]
 *     [ $projection_on_f0_f1_axis, $chi_square, $deg_of_freedom, $n_codon, $id ]
 *
 *------------------------------------------------------------------------------
 *
 *  For a given amino acid, the relative frequency of codon i at point x is
 *  f(i,x).
 *
 *  Define frequencies f0(i) and f1(i) that we want to match at x = 0 and 1.
 *  Typically these will be the modal and high expression usage frequencies.
 *  Define the relative weight of codon i at point x as:
 *
 *     w(i,x) = f0(i) * exp( k(i) * x )
 *
 *  where f0(i) is the frequency of codon i at x = 0.  To get the frequency
 *  of codon i at point x we normalize by the sum of the weights at x:
 *
 *     f(i,x) = w(i,x) / sum_over_j( w(j,x) )
 *
 *  At x = 0,
 *
 *     f(i,0) = w(i,0) / sum_over_i( w(i,0) )
 *            = f0(i) * exp( k(i) * 0 ) / sum_over_j( f0(j) * exp( k(j) * 0 ) )
 *            = f0(i) * exp( 0 ) / sum_over_j( f0(j) * exp( 0 ) )
 *            = f0(i) * 1 / sum_over_j( f0(j) * 1 )
 *            = f0(i) / sum_over_j( f0(j) )
 *            = f0(i) / 1
 *            = f0(i)
 *
 *  as it should.  Note that k(i) does not enter.  At x = 1,
 *
 *     f(i,1)   f0(i) * exp( k(i) )
 *     ------ = -------------------
 *     f(j,1)   f0(j) * exp( k(j) )
 *
 *  or
 *
 *     exp( k(i) )   f(i,1) / f0(i)   f1(i) / f0(i)
 *     ----------- = -------------- = -------------
 *     exp( k(j) )   f(j,1) / f0(j)   f1(j) / f0(j)
 *
 *  where f1(i) = f(i,1).
 *
 *  The values of k are underdetermined within a constant multiplier, but this
 *  form suggests that it is natural to define:
 *
 *     exp( k(i) ) = f1(i) / f0(i)
 *
 *  or
 *
 *     k(i) = ln( f1(i) / f0(i) )
 *
 *-----------------------------------------------------------------------------
 *
 *  sub project_on_freq_vector_by_chi_sqr_2
 *  {
 *      my $f0 = shift;                               #  Frequencies at x = 0
 *      my $f1 = shift;                               #  Frequencies at x = 1
 *      my $opts = ref $_[0] eq 'HASH' ? shift : {};  #  Options
 *      return () if ! @_;
 *  
 *      return undef if zero_vector( subtract_points( $f1, $f0 ) );
 *  
 *      #  Compute exponential coeficient k for each codon for each amino acid:
 *  
 *      my $k = k_from_f0_and_f1( $f0, $f1 );
 *  
 *      #  For each set of counts in the arguements, find the position along the
 *      #  $freq0 to $freq1 vector that has smallest chi square.  Figuring out
 *      #  the nature of the input list is a bit convoluted.
 *      #
 *      #  ! ref $_[0]->[0]->[0]              -->  $_[0] is   $counts
 *      #  @{$_[0]} == 2 && ! ref $_[0]->[1]  -->  $_[0] is [ $counts, $id ]
 *
 *      my @projections = map { project_by_min_chi_sqr_2( $f0, $k, $_ ) }
 *                        ( ! ref $_[0]->[0]->[0]              ? @_
 *                        :  @{$_[0]} == 2 && ! ref $_[0]->[1] ? @_
 *                        :                                      @{$_[0]}
 *                        );
 *  
 *      wantarray ? @projections : \@projections
 *  }
 *
 *
 *------------------------------------------------------------------------------
 *  If $counts includes an id, i.e., $counts = [ [ [ ... ], ... ], $id ]:
 *
 *  [ $x, $chisqr, $df, $ncodon, $id ] = project_by_min_chi_sqr_2( $f0, $k, $counts )
 *
 *  If $counts does not include an id, i.e., $counts = [ [ ... ], ... ]:
 *
 *  [ $x, $chisqr, $df, $ncodon ]      = project_by_min_chi_sqr_2( $f0, $k, $counts )
 *------------------------------------------------------------------------------
 */

x_pval_t
project_by_min_chisqr( codon_freq_t f0, k_vector_t k, codon_cnt_t counts, int maxlen )
{
    #define MIN_X  -5.0
    #define MAX_X   5.0
    #define N_STEP 20
    #define MIN_DX 0.0005

    x_pval_t    result;
    chisqr_n_t  best = { 1e99, 0, 0 };
    chisqr_n_t  chisqr;
    double      x, dx;
    double      best_x = 0;
    int         i;

    /*  Grid search centered on 0  */

    dx = ( MAX_X - MIN_X ) / N_STEP;
    x  =   MIN_X + dx;
    for ( i = 1; i < N_STEP; i++ )
    {
        chisqr = chisqr_at_x( counts, f0, k, x );
        if ( chisqr.chisqr < best.chisqr )
        {
            best_x      = x;
            best.chisqr = chisqr.chisqr;
            best.df     = chisqr.df;
            best.n      = chisqr.n;
        }
        x += dx;
    }

    /*  Divide and conquer search centered on current best point  */

    while ( dx > MIN_DX )
    {
        dx *= 0.5;
        x = best_x - dx;
        chisqr = chisqr_at_x( counts, f0, k, x );
        if ( chisqr.chisqr < best.chisqr )
        {
            best_x      = x;
            best.chisqr = chisqr.chisqr;
            best.df     = chisqr.df;
            best.n      = chisqr.n;
        }
        x = best_x + dx;
        chisqr = chisqr_at_x( counts, f0, k, x );
        if ( chisqr.chisqr < best.chisqr )
        {
            best_x      = x;
            best.chisqr = chisqr.chisqr;
            best.df     = chisqr.df;
            best.n      = chisqr.n;
        }
    }

    /* return best; */

    result.x = best_x;
    result.p = chisqr_prob_w_maxn( best.chisqr, best.df, best.n, maxlen );

    return result;
}


/*-----------------------------------------------------------------------------
 *  For a given point x along the (extrapolated) line from f0 to f1, find
 *  the codon frequencies.  The exponential coeficients k can be computed
 *  from f0 and f1 with k_from_f0_and_f1( f0, f1 ).
 *-----------------------------------------------------------------------------
 */

codon_freq_t
freqs_at_x( codon_freq_t f0, k_vector_t k, double x )
{
    codon_freq_t f;
    double *f0aa;
    double *kaa;
    double *fraa;
    double fr;
    double sum_fr;
    double sum_fr_2;
    double min_f = 0.0001;       /* No frequencies less than min_f  */
    int aa;
    int offset;
    int i;
    int ilim;
    int renorm;

    /*  Compute freqs for each amino acid and codon:  */

    for ( aa = 0; aa <= 17; aa++ )
    {
        offset   = aa_offset[ aa ];
        f0aa     = &f0.freq[ offset ];
        kaa      = &k.k[ offset ];
        fraa     = &f.freq[ offset ];
        sum_fr   = 0;
        sum_fr_2 = 0;
        ilim     = aa_n_cdn[aa];
        renorm   = 0;

        /*  Calculate unnormalized frequencies for the codons of the amino
         *  acid.
         */
        for ( i = 0; i < ilim; i++ )
        {
            sum_fr += *(fraa + i) = *(f0aa + i) * exp( *(kaa + i) * x );
        }

        /*  Normalize frequencies, recording minimum found  */

        for ( i = 0; i < ilim; i++ )
        {
            fr = *(fraa + i) /= sum_fr;
            if ( fr < min_f ) { fr = *(fraa + i) = min_f; renorm = 1; }
            sum_fr_2 += fr;
        }

        /*  If any frequencies were below min_f, normalize again:  */

        if ( renorm )
        {
            for ( i = 0; i < ilim; i++ ) { *(fraa + i) /= sum_fr_2; }
        }

    }

    return f;
}


/*-----------------------------------------------------------------------------
 *  Compute the exponential coeficient k for each codon for each amino acid:
 *  (does not check than frequencies are greater than 0)
 *-----------------------------------------------------------------------------
 */

k_vector_t
k_from_f0_and_f1( codon_freq_t f0, codon_freq_t f1 )
{
    k_vector_t k;
    int i;
    int ilim = aa_offset[ 18 ]; /* first unused codon */
    for ( i = 0; i <= ilim; i++ ) { k.k[i] = log( f1.freq[i] / f0.freq[i] ); }

    k.naa = 18;
    return k;
}


/*-----------------------------------------------------------------------------
 *  P-value for a gene codon usage for a given x value.
 *-----------------------------------------------------------------------------
 */
double
pval_at_x( codon_cnt_t counts, int maxlen, codon_freq_t f0, k_vector_t k, double x )
{
    chisqr_n_t  chi2_n;
    chi2_n = count_vs_freq_chisqr( counts, freqs_at_x( f0, k, x ) );
    return chisqr_prob_w_maxn( chi2_n.chisqr, chi2_n.df, chi2_n.n, maxlen );
}


/*-----------------------------------------------------------------------------
 *  chi-square for a gene codon usage for a given x value.
 *-----------------------------------------------------------------------------
 */
chisqr_n_t
chisqr_at_x( codon_cnt_t counts, codon_freq_t f0, k_vector_t k, double x )
{
    return count_vs_freq_chisqr( counts, freqs_at_x( f0, k, x ) );
}


/*-----------------------------------------------------------------------------
 *  For a given point x along the (extrapolated) line from f0 to f1, find
 *  the codon frequencies.  (When used multiple times for the same f0 and
 *  f1, it is more efficient to save the exponentical coeficients k.)
 *-----------------------------------------------------------------------------
 */

codon_freq_t
freqs_at_x_given_f0_and_f1( codon_freq_t f0, codon_freq_t f1, double x )
{
    return freqs_at_x( f0, k_from_f0_and_f1( f0, f1 ), x );
}


/*=============================================================================
 *  Read and write counts and frequencies
 *-----------------------------------------------------------------------------
 *  Read a line of codon counts from a file.
 *
 *     codon_cnt_t read_codon_counts( FILE * fp )
 *
 *  Counts are space separated. If there is a trailing id, it is tab separated.
 *  Test for success by result.naa >= 18.
 *-----------------------------------------------------------------------------
 */
 
codon_cnt_t
read_codon_counts( FILE * fp )
{
    codon_cnt_t  result;
    char         line[READBUF];
    char        *read;
    char        *labelptr = (char *) NULL;

    read = fgets( line, READBUF, fp );
    if ( read ) { result = parse_codon_counts( line, &labelptr ); }
    else        { result.naa = 0; }

    return result;
}


/*-----------------------------------------------------------------------------
 *  Read a line of codon counts from a file, setting a pointer to the label.
 *
 *     codon_cnt_t read_codon_counts_2( FILE * fp, char * label, int max_lbl_len )
 *
 *  Counts are space separated. If there is a trailing id, it is tab separated.
 *  Test for success by result.naa >= 18.
 *-----------------------------------------------------------------------------
 */
 
codon_cnt_t
read_codon_counts_2( FILE * fp, char * label, int max_lbl_len )
{
    codon_cnt_t  result;
    char         line[READBUF];
    char        *read;
    char        *labelptr = (char *) NULL;

    read = fgets( line, READBUF, fp );
    if ( read )
    {
        result = parse_codon_counts( line, &labelptr );
        if ( labelptr && *labelptr )
        {
            scopyn( labelptr, label, max_lbl_len );
        }
        else
        {
            *label = '\0';
        }
    }
    else
    {
        result.naa = 0;
        *label = '\0';
    }

    return result;
}


/*-----------------------------------------------------------------------------
 *  A simple string copy for the label.  Adds "\n\0" if it must truncate.
 *
 *     scopyn( char * s1, char * s2, int max )
 *
 *-----------------------------------------------------------------------------
 */
void
scopyn( char * s1, char * s2, int max )
{
    if ( max > 0 )
    {
        while ( *s2++ = *s1++ )
        {
            if ( ! --max ) { *--s2 = '\0'; *--s2 = '\n'; break; }
        }
    }
}

/*-----------------------------------------------------------------------------
 *  Read codon counts from a string.
 *
 *     codon_cnt_t parse_codon_counts( char * line, char ** labelptr )
 *
 *  Counts are space separated. If there is a trailing id, it is tab separated.
 *  Test for success by result.naa >= 18.
 *-----------------------------------------------------------------------------
 */
codon_cnt_t
parse_codon_counts( char * line, char ** labelptr )
{
    codon_cnt_t  result;
    char        *cp = line;
    int         *cnt;
    int          n, aa, cdn;

    result.naa = 0;
    cnt = &result.cnt[0];

    for ( aa = 0; aa < 20; aa++ )
    {
        for ( cdn = 0; cdn < aa_n_cdn[aa]; cdn++ )
        {
            n = s_next_count( &cp );
            if ( n < 0 )
            {
                *labelptr = ( n == -3 ) ? cp : (char *) NULL;
                return result;
            }
            *cnt++ = n;
        }
        result.naa = aa + 1;     /* finished the amino acid */
    }

    *labelptr = ( *cp == '\t' ) ? cp+1 : (char *) NULL;
    return result;
}


/* Return values:
 *  >= 0 count
 *    -1 EOF
 *    -2 newline   # consumed
 *    -3 tab       # consumed
 *    -4 nondigit  # pushed back
 */
int
s_next_count( char **cpp )
{
    int c, i;
    c = s_flush_count_space( cpp );
    if ( c <=  0 )        { return -1; }
    if ( c == '\n' )      { return -2; }
    if ( c == '\t' )      { return -3; }
    if ( ! isdigit( c ) ) { return -4; }

    i = 0;
    while ( ( c = s_get_digit( cpp ) ) > 0 ) { i = 10 * i + ( c - '0' ); }
    return i;
}


int
s_get_digit( char ** cpp )
{
    int c;
    c = **cpp;
    if ( ! isdigit( c ) ) { return -1; }
    (*cpp)++;   /* digits are consumed */
    return c;
}


/*  Consume blank space, returning the first non-space character, tab or
 *  newline. Tab and new line are consumed. Other terminating characters
 *  are pushed back on the input.
 */
int
s_flush_count_space( char **cpp )
{
    int c;
    while ( ( ( c = **cpp ) > 0 ) && isspace( c ) )
    {
        (*cpp)++;
        if ( c == '\n' || c == '\t' ) { break; }  /* spaces that end scan */
    }
    return c;
}


/*-----------------------------------------------------------------------------
 *  Read a line of codon frequencies from a file.
 *
 *     codon_freq_t read_codon_freqs( FILE * fp )
 *
 *  Frequencies are separated by '|' or ','. If there is a leading score, it
 *  is tab separated.  If there is a trailing id, it is tab separated.
 *  Test for success by result.naa >= 18.
 *-----------------------------------------------------------------------------
 */
codon_freq_t
read_codon_freqs( FILE * fp )
{
    codon_freq_t  result;
    char          line[READBUF];
    char         *read;

    read = fgets( line, READBUF, fp );
    if ( read ) { result = parse_codon_freqs( line ); }
    else        { result.naa = 0; }

    return result;
}


/*-----------------------------------------------------------------------------
 *  Read codon frequencies from a string.
 *
 *     codon_freq_t parse_codon_freqs( char * line )
 *
 *  Frequencies are separated by '|' or ','. If there is a leading score, it
 *  is tab separated.  If there is a trailing id, it is tab separated.
 *  Test for success by result.naa >= 18.
 *-----------------------------------------------------------------------------
 */
codon_freq_t
parse_codon_freqs( char * line )
{
    codon_freq_t  result;
    char         *cp = line;
    double       *freq, f;
    int           c, aa, cdn;

    /* Figure out if there is a leading score to be skipped.  */

    while ( ( c = *cp ) && c == ' ' ) { cp++; }   /* leading blanks */
    while ( ( c = *cp ) && ( isdigit( c ) || c == '.' || c == '-' || c == 'e' ) ) { cp++; } /* number */
    while ( ( c = *cp ) && c == ' ' ) { cp++; }   /* trailing blanks */

    if ( ( c = *cp ) && c == '\t' ) { cp++; }       /* score found and removed */
    else                            { cp = line; }  /* no score; reset pointer */

    result.naa = 0;
    freq = &result.freq[0];
    for ( aa = 0; aa < 20; aa++ )
    {
        for ( cdn = 0; cdn < aa_n_cdn[aa]; cdn++ )
        {
            f = s_next_frequency( &cp );
            if ( f < 0 ) { return result; }
            *freq++ = f;
        }
        result.naa = aa + 1;     /* finished the amino acid */
    }

    return result;
}


/* Return values:
 *  >= 0 frequency
 *    -1 EOF
 *    -2 newline   # consumed
 *    -3 tab       # consumed
 *    -4 nondigit  # pushed back
 */
double
s_next_frequency( char **cpp )
{
    double  x;
    int     c;
    c = s_flush_freq_div( cpp );
    if ( c <=  0 )        { return -1; }
    if ( c == '\n' )      { return -2; }
    if ( c == '\t' )      { return -3; }
    if ( ! isdigit( c ) ) { return -4; }

    x = 0;
    while ( ( c = s_get_digit( cpp ) ) > 0 ) { x = 10 * x + ( c - '0' ); }
    if ( **cpp == '.' )
    {
        double y = 0.1;
        (*cpp)++;  /* consume the decimal point */
        while ( ( c = s_get_digit( cpp ) ) > 0 )
        {
            x += y * ( c - '0' );
            y *= 0.1;
        }
    }
    return x;
}


/*  This will silently consume multiple dividers.  Perhaps they should
 *  be contrived to be zero values, but for now ....
 */
int
s_flush_freq_div( char ** cpp )
{
    int c;
    while ( ( ( c = **cpp ) > 0 )
         && ( isspace( c ) || c == '|' || c == ',' )
          )
    {
        (*cpp)++;
        if ( c == '\n' || c == '\t' ) { break; }  /* spaces that end scan */
    }
    return c;
}


/*-----------------------------------------------------------------------------
 *  Print a line of codon counts to a file.
 *
 *     void  fprint_codon_counts( FILE * fp, codon_cnt_t counts )
 *
 *-----------------------------------------------------------------------------
 */
void
fprint_codon_counts( FILE * fp, codon_cnt_t counts )
{
    int  aa;
    int  cdn;
    int  offset = 0;

    for ( aa = 0; aa < counts.naa; aa++ )
    {
        if ( aa > 0 ) { putc( ' ', fp ); }
        for ( cdn = 0; cdn < aa_n_cdn[aa]; cdn++ )
        {
            if ( aa > 0 || cdn > 0 ) { putc( ' ', fp ); }
            fprintf( fp, "%d", counts.cnt[offset++] );
        }
    }
    putc( '\n', fp );
}


/*-----------------------------------------------------------------------------
 *  Print a line of codon frequencies to a file.
 *
 *     void  fprint_codon_freqs( FILE * fp, codon_freq_t freqs )
 *
 *-----------------------------------------------------------------------------
 */
void
fprint_codon_freqs( FILE * fp, codon_freq_t freqs )
{
    int  aa;
    int  cdn;
    int  offset = 0;

    for ( aa = 0; aa < freqs.naa; aa++ )
    {
        if ( aa > 0 ) { putc( '|', fp ); }
        for ( cdn = 0; cdn < aa_n_cdn[aa]; cdn++ )
        {
            if ( cdn > 0 ) { putc( ',', fp ); }
            fprintf( fp, "%.5f", freqs.freq[offset++] );
        }
    }
    putc( '\n', fp );
}


/*------------------------------------------------------------------------------
 *  Redefine the number of codons for each amino acid.  The expected order is:
 *
 *     A  C  D  E  F  G  H  I  K  L  N  P  Q  R  S  T  V  Y  M  W
 *
 *  Values for M and W are optional.  Values are saved in the global variables:
 *
 *     int aa_n_cdn[20];
 *     int aa_offset[20];
 *
 *-----------------------------------------------------------------------------
 */
int
set_codons_per_aa( int naa, int * new_aa_n_cdn )
{
    int aa;
    int offset = 0;
    int okay = 1;

    if ( naa < 18 || naa > 20 )
    {
        fprintf( stderr, "Warning: set_codons_per_aa called with %i amino acids.\n", naa );
        if ( naa > 20 ) { naa = 20; }
        okay = 0;
    }

    for ( aa = 0; aa < naa; aa++ )
    {
        aa_offset[ aa ] = offset;
        offset += aa_n_cdn[ aa ] = new_aa_n_cdn[ aa ];
    }

    for ( ; aa < 20; aa++ )
    {
        aa_offset[ aa ] = offset++;
        aa_n_cdn[ aa ] = 1;
    }

    if ( offset > 64 )
    {
        fprintf( stderr, "Warning: set_codons_per_aa called with %i codons.\n", offset );
        okay = 0;
    }

    return okay;
}


/*  Set codons for each amino acid to the defaults.  */
int
reset_codons_per_aa()
{
    /*                        A  C  D  E  F  G  H  I  K  L  N  P  Q  R  S  T  V  Y  M  W */
    int  new_aa_n_cdn[20] = { 4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 2, 4, 2, 6, 6, 4, 4, 2, 1, 1 };
    return set_codons_per_aa( 20, new_aa_n_cdn );
}


/*=============================================================================
 *  Statistical routines
 *=============================================================================
 *  chi-square value for n counts matching a set of n frequencies.
 *-----------------------------------------------------------------------------
 */
chisqr_n_t
chisqr( double * freq, int * cnt, int n )
{
    chisqr_n_t  result = { 0, -1, 0 };
    double      sum = 0;
    double      ttlcnt_d, diff, ex, fr;
    int         ttlcnt = 0;
    int         i;
    int         df = -1;
    int        *c  = cnt;

    for ( i = 0; i < n; i++ ) { ttlcnt += *c++; }
    if ( ttlcnt <= 0 ) { return result; }  /* no counts */

    ttlcnt_d = (double) ttlcnt;
    for ( i = 0; i < n; i++ )
    {
        fr = *freq++;
        if ( fr <= 0 ) { continue; } /* should diagnose error if counts */
        ex = fr * ttlcnt_d;
        diff = *cnt++ - ex;
        sum += diff * diff / ex;
        df++;
    }

    result.chisqr = sum;
    result.df     = df;
    result.n      = ttlcnt;

    return result;
}


/*------------------------------------------------------------------------------
 *  A chi square test with a limit on the total observed. If the total
 *  observations exceeds maxn, then the chi square value is scaled down
 *  before the p-value calculation.
 *-----------------------------------------------------------------------------
 */
double
chisqr_prob_w_maxn( double chi2, int df, int n, int maxn )
{
    if ( ( df < 1 ) || ( n < 1 ) ) { return 1; }
    if ( maxn && n > maxn ) { chi2 *= (double) maxn / (double) n; }
    return chisqr_prob( chi2, df );
}


/*-----------------------------------------------------------------------------
 *  Probability of a chi square value greater than or equal to chisqr
 *
 *  Based on:
 *
 *  Zelen, M. and Severo, N. C. (1965).  Probability functions.  In
 *  Handbook of Mathematical Functions, Abramowitz, M. and Stegun,
 *  I. A., eds. (New York: Dover Publications), pp. 925-995.
 *
 *  Programmed in C by Gary Olsen
 *
 *     double p_value = chisqr_prob( double chisqr, int df )
 *
 *-----------------------------------------------------------------------------
 */

double
chisqr_prob( double chisqr, int df )
{
    double sum, delta, denom1, denom2, prob_not;
    int    i;

    if ( ( chisqr <  0 ) || ( df <  0 ) ) { return -1; }
    if ( ( chisqr == 0 ) || ( df == 0 ) ) { return  1; }
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


