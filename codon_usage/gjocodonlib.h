/* gjocodonlib.h */

#define READBUF 256000

/*  Much of the code relies on "packaged codon counts" and "packaged
 *  codon frequencies. The codons per amino acid can be changed with
 *  set_codons_per_aa() and reset_codons_per_aa().
 *
 *                     A  C  D  E  F  G  H  I  K  L  N  P  Q  R  S  T  V  Y  M  W */
int  aa_n_cdn[20]  = { 4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 2, 4, 2, 6, 6, 4, 4, 2, 1, 1 };
int  aa_offset[20] = { 0, 4, 6, 8,10,12,16,18,21,23,29,31,35,37,43,49,53,57,59,60 };

char one_ltr_aa[20] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                        'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'M', 'W'
                      };

typedef struct { int naa; int val[64]; }     codon_int_t;
typedef struct { int naa; int cnt[64]; }     codon_cnt_t;

typedef struct { int naa; double val[64];  } codon_dbl_t;
typedef struct { int naa; double freq[64]; } codon_freq_t;

typedef struct { int naa; double k[64];    } k_vector_t;

/*  The codon order in the perl programs in gjocodonlib.pm:
 *
 *  %amino_acid_codons_DNA = (
 *           A  => [ qw( GCA GCG GCT GCC ) ],
 *           C  => [ qw( TGT TGC ) ],
 *           D  => [ qw( GAT GAC ) ],
 *           E  => [ qw( GAA GAG ) ],
 *           F  => [ qw( TTT TTC ) ],
 *           G  => [ qw( GGA GGG GGT GGC ) ],
 *           H  => [ qw( CAT CAC ) ],
 *           I  => [ qw( ATA ATT ATC ) ],
 *           K  => [ qw( AAA AAG ) ],
 *           L  => [ qw( TTA TTG CTA CTG CTT CTC ) ],
 *           N  => [ qw( AAT AAC ) ],
 *           P  => [ qw( CCA CCG CCT CCC ) ],
 *           Q  => [ qw( CAA CAG ) ],
 *           R  => [ qw( AGA AGG CGA CGG CGT CGC ) ],
 *           S  => [ qw( AGT AGC TCA TCG TCT TCC ) ],
 *           T  => [ qw( ACA ACG ACT ACC ) ],
 *           V  => [ qw( GTA GTG GTT GTC ) ],
 *           Y  => [ qw( TAT TAC ) ],
 *           M  => [ qw( ATG ) ],
 *           W  => [ qw( TGG ) ]
 *         );
 */

/* some probability and chi-square data structures */

typedef struct { double x; double p; }                      x_pval_t;
typedef struct { double x; double chisqr; int df; int n; }  x_chisqr_t;
typedef struct { double chisqr; int df; }                   chisqr_t;
typedef struct { double chisqr; int df; int n; }            chisqr_n_t;

/* function prototypes */

double        count_vs_freq_pval( codon_cnt_t cnt, codon_freq_t freq, int maxlen );
chisqr_n_t    count_vs_freq_chisqr( codon_cnt_t cnt, codon_freq_t freq );

x_pval_t      project_by_min_chisqr( codon_freq_t f0, k_vector_t k, codon_cnt_t counts, int maxlen );
codon_freq_t  freqs_at_x( codon_freq_t f0, k_vector_t k, double x );
k_vector_t    k_from_f0_and_f1( codon_freq_t f0, codon_freq_t f1 );
double        pval_at_x( codon_cnt_t counts, int maxlen, codon_freq_t f0, k_vector_t k, double x );
chisqr_n_t    chisqr_at_x( codon_cnt_t counts, codon_freq_t f0, k_vector_t k, double x );
codon_freq_t  freqs_at_x_given_f0_and_f1( codon_freq_t f0, codon_freq_t f1, double x );

codon_cnt_t   read_codon_counts( FILE * fp );
codon_cnt_t   read_codon_counts_2( FILE * fp, char * label, int max_lbl_len );
codon_cnt_t   parse_codon_counts( char * cp, char ** labelp );
int           s_next_count( char ** cpp );
int           s_get_digit( char ** cpp );
int           s_flush_count_space( char ** cpp );
void          scopyn( char * s1, char * s2, int max );

codon_freq_t  read_codon_freqs( FILE * fp );
codon_freq_t  parse_codon_freqs( char * cp );
double        s_next_frequency( char ** cpp );
int           s_flush_freq_div( char ** cpp );

void          fprint_codon_counts( FILE * fp, codon_cnt_t counts );
void          fprint_codon_freqs( FILE * fp, codon_freq_t freqs );

int           set_codons_per_aa( int naa, int * new_aa_n_cdn );
int           reset_codons_per_aa();

chisqr_n_t    chisqr( double * freq, int * cnt, int n );
double        chisqr_prob( double chisqr, int df );
double        chisqr_prob_w_maxn( double chi2, int df, int n, int maxn );
