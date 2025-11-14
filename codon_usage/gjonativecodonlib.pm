package gjonativecodonlib;
#
#  A library of functions for finding native (modal and high expression)
#  codon usage of a genome.
#
use strict;
use bidir_best_hits qw( bbh );
use gjocodonlib     qw( codon_freq_distance
                        count_vs_freq_chi_sqr
                        seq_codon_count_package
                        report_frequencies
                        modal_codon_usage
                        project_on_freq_vector_by_chi_sqr_2
                        report_counts
                      );
use gjoseqlib       qw( translate_seq );
use gjostat         qw( chisqr_prob mean_stddev );
use IPC::Open2      qw( open2 );
# use Data::Dumper;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( native_codon_usage );

my $usage = <<'End_of_Comments';
Find the modal and high expression codon usage of a genome.

   @mode_label_pairs = native_codon_usage( \%params )

The first mode returned is that of the genome.  By default, the second mode
returned is the final estimate.  If all_he_modes is true, then the series of
estimates is included.

Required parameter:

    dna => \@dna_entries  # required

Optional parameters:

    all_he_modes     => bool           # D = 0       # show successive high expr approximations
    average          => bool           # D = undef   # include the average codon usage
    bbh_blast_opts   => opts           # D = ''      # none
    bbh_coverage     => frac           # D = 0.70    # bbh covers 70% of gene
    bbh_e_value      => e_value        # D = 1e-5    # bbh match at least 1e-5 (could relax a little)
    bbh_identity     => frac           # D = 0.20    # bbh 20% amino acid identity
    bbh_positives    => frac           # D = 0.30    # bbh 30% amino acids with positive score
    bias_stats       => bool           # D = undef   # don't do bias statistics
    counts_file      => file           # D = ''      # none
    genome_title     => string         # D = ''      # prefix to each mode label
    log              => file_handle    # D = STDERR  # (if log_label is supplied)
    log_label        =>                # D = ''      # no log
    he_ref_module    => perl_module    # D ='high_expression_ref_ab';
    match_p_value    => p_value        # D = 0.10;   # a match it mode is P > 0.1
    max_he_decline   => fraction       # D = 0.20    # < 20% decline in finding original he candidates
    max_iter         => int            # D = 2       # original estimate and 2 refinements
    max_keep         => int or frac    # D = 0.10    # keep up to 10% of genome
    max_mode_overlap => frac           # D = 0.80    # < 80% of genes matching he can match mode
    min_he_genes     => int            # D = 20      # don't estimate with < 20 bbh genes
    min_nonmodal     => frac           # D = 0.20    # 20% of candidates differ from genome mode
    min_x_value      => frac           # D = 0.50    # >%50 of way from mode to he estimate
    mode_exponent    => real           # D = 0.30    # exponent of P-value in mode calculation
    omit_p_value     => p_value        # D = 0.10    # omit from high exp calc if P > 0.1 to mode
    p_val_max_len    => int            # D = undef;  # no limit on length in chi-square
    reach_p_value    => p_value        # D = 0.05    # allow direction to drift

Special parameters that change the function, or return values:

   algorithm => bool  # if true, the function returns the algorthm description
   usage     => bool  # if true, the function returns this usage description
   version   => bool  # if true, the function returns the version string

If none of these are true, then the analysis is performed, and the following
assignments are made:

   $param->{ algorithm } = \$algorithm_description
   $param->{ usage }     = \$usage_description
   $param->{ version }   =  $version_string

End_of_Comments

my $algorithm = <<'End_of_Algorithm';
Details of the Algorithm and Parameters:

Find the native (modal and high expression) codon usages of a genome.
The algorithm can be summarized as:

1.  Find the modal codon usage of the genome.

    Relevant paramter:

        mode_exponent -
             Exponent of P-value in optimizing the mode.

2.  Identify a set of candidate high expression genes by finding the genes
    that are bidirectional best hits to presumptively abundant proteins from
    a reference genome.

    Relevant paramters:

        bbh_blast_opts -
             Blast program options (note: most blast options are ignored).
                Perhaps the most useful is '-a 4' for using 4 threads.
        bbh_coverage -
             Minimum fraction of query and subject sequence lengths covered
                in blast matches
        bbh_e_value -
             Maximum blast E-value
        bbh_identity -
             Minimum (fraction) sequence identity in blast matches
        bbh_positives -
             Minimun (fraction) amino acid matches with positive score in
                blast search
        he_ref_module -
             Perl module with reference genome and high expression gene list

3.  Use a set of candidate high expression genes (from step 2 or step 5) to
    estimate the high expression codon usage of the genome. In detail:

    3a. Remove from the candidate high expression genes those that are too
        similar to the overall mode (from step 1).

        Relevant paramters:

            min_nonmodal -
                 Fraction of candidate high expression genes that must differ
                    significantly from the genome mode for fillering to occur
            omit_p_value -
                 Genes that match mode with P > omit_p_value are usually
                    omitted from high expression estimate.  Thus,l
                    omit_p_value => 1 wil use all of the genes.
            p_val_max_len -
                 Length limit in calculating match of protein to the modal
                    codon usage. Longer proteins have their chi-square scaled
                    down before calculating the P-value, so large genes appear
                    to match better.

    3b. Find the modal usage of the remaining candidate high expression genes.

        Relevant paramter:

            min_he_genes -
                 Number of candidate high expression genes for the mode to
                    be computed

4.  Terminate if a stopping condition is reached.

        Relevant paramters:

            max_iter -
                 Number of iterations of refinement, that is cycles of steps
                    5, 3, and 4
            min_nonmodal -
                 Fraction of candidate high expression genes that must be
                    distinct (significantly different) from the modal codon
                    usage for the iterations to continue

5.  Otherwise, produce a new set of candidate high expression genes. In detail:

    5a. Limit to genes that are within "reach P-value" of current
        "mode to high expression" axis

        Relevant paramters:

            p_val_max_len -
                 Length limit in calculating match of protein to the modal
                    codon usage. Longer proteins have their chi-square scaled
                    down before calculating the P-value, so large genes appear
                    to match better.
            reach_p_value -
                 Minimum P-value of a match to the mode to high expression
                    axis for the gene to be included in the next estimate.
                    This defines which genes are within "reach" for influencing
                    the next mode calculation. The goal is to allow migration
                    of the estimate within a constrained space.

    5b. Limit to genes with x value greater than or equal to "min x value".

        Relevant paramter:

            min_x_value -
                 Minium value of x for inclusion as candidate high expression
                    gene. For any given iteration, either max_keep or
                    min_x_value will be limiting, as they both monitor the
                    x value of the genes.

    5c. Limit to "max keep" highest x values (most high expression like).

        Relevant paramter:

            max_keep -
                 Maxium number of candidate high expression genes kept for the
                    next estimate (before removing those that match the mode).
                    For any given iteration, either max_keep or min_x_value
                    will be limiting, as they both monitor the x value of the
                    genes.

    5d. Go to step 3. This means that max keep (max_keep value) is enforced
        before removing the mode-matching genes.

End_of_Algorithm

#===============================================================================
#  Find the genome modal codon usage, and one or more estimates of the high
#  expression codon usage.
#
#  ( $genome_mode, @high_exp_mode_ests ) = native_codon_usage( \%options )
#
#  Each mode is a mode-label pair.  With valid data, the genome mode should
#  always be returned.  The caller is responsible for obvious filtering,
#  such as minimum DNA sequence length.
#===============================================================================

sub native_codon_usage
{
    $_[0] && ref $_[0] eq 'HASH' or return undef;

    my $version = '1.00';
    my $param = $_[0];

    #  Special requests:

    return $version   if ( $param->{ version }   && $param->{ version }   ne  $version );
    return $usage     if ( $param->{ usage }     && $param->{ usage }     ne \$usage );
    return $algorithm if ( $param->{ algorithm } && $param->{ algorithm } ne \$algorithm );

    #  Otherwise, link the information into the calling paramter array:

    $param->{ version }   =  $version;    #  Return version version
    $param->{ usage }     = \$usage;      #  Return usage information
    $param->{ algorithm } = \$algorithm;  #  Return algorithm information

    # use high_expression_bacteria qw( &ref_aa &ref_he )   # embedded in the code

    #  \$(\S+)(\s*)=
    #  $\1\2= defined( $param->{ \1 } )\2? $param->{ \1 }\2:

    #  Parameters that the user can change:

    my $all_he_modes     = defined( $param->{ all_he_modes } )     ? $param->{ all_he_modes }     :     0;      # include successive approximations
    my $average          = defined( $param->{ average } )          ? $param->{ average }          :     0;      # include average as first frequency
    my $bbh_blast_opts   = defined( $param->{ bbh_blast_opts } )   ? $param->{ bbh_blast_opts }   :    '';
    my $bbh_coverage     = defined( $param->{ bbh_coverage } )     ? $param->{ bbh_coverage }     :     0.70;   # bbh covers 70% of gene
    my $bbh_e_value      = defined( $param->{ bbh_e_value } )      ? $param->{ bbh_e_value }      :  1e-5;      # bbh match at least 1e-5 (could relax a little)
    my $bbh_identity     = defined( $param->{ bbh_identity } )     ? $param->{ bbh_identity }     :     0.20;   # bbh 20% amino acid identity
    my $bbh_positives    = defined( $param->{ bbh_positives } )    ? $param->{ bbh_positives }    :     0.30;   # bbh 30% amino acids with positive score
    my $bias_stats       = defined( $param->{ bias_stats } )       ? $param->{ bias_stats }       : undef;      # don't do bias statistics
    my $counts_file      = defined( $param->{ counts_file } )      ? $param->{ counts_file }      :    '';      # don't save the counts
    my $dna              = defined( $param->{ dna } )              ? $param->{ dna }              :    [];
    my $genome_title     = defined( $param->{ genome_title } )     ? $param->{ genome_title }     :    '';      # prefix for codon frequency labels
    my $he_ref_module    = defined( $param->{ he_ref_module } )    ? $param->{ he_ref_module }    : 'high_expression_ref_ab';
    my $initial_he_ids   = defined( $param->{ initial_he_ids } )   ? $param->{ initial_he_ids }   : undef;
    my $log              = defined( $param->{ log } )              ? $param->{ log }              :    '';      # file handle for log
    my $log_label        = defined( $param->{ log_label } )        ? $param->{ log_label }        :    '';
    my $match_p_value    = defined( $param->{ match_p_value } )    ? $param->{ match_p_value }    :     0.10;   #
    my $max_he_decline   = defined( $param->{ max_he_decline } )   ? $param->{ max_he_decline }   :     0.20;   # < 20% decline in finding original he candidates
    my $max_iter         = defined( $param->{ max_iter } )         ? $param->{ max_iter }         :     2;      # original estimate and 2 refinements
    my $max_keep         = defined( $param->{ max_keep } )         ? $param->{ max_keep }         :     0.10;   # 10% of the genome
    my $max_mode_overlap = defined( $param->{ max_mode_overlap } ) ? $param->{ max_mode_overlap } :     0.80;   # < 80% of genes matching he can match mode
    my $min_he_genes     = defined( $param->{ min_he_genes } )     ? $param->{ min_he_genes }     :    20;      # don't estimate with < 20 bbh genes
    my $min_nonmodal     = defined( $param->{ min_nonmodal } )     ? $param->{ min_nonmodal }     :     0.20;   # 20% of candidates differ from genome mode
    my $min_x_value      = defined( $param->{ min_x_value } )      ? $param->{ min_x_value }      :     0.50;   #
    my $mode_exponent    = defined( $param->{ mode_exponent } )    ? $param->{ mode_exponent }    :     0.30;   # exponent of P-value in mode calculation
    my $omit_p_value     = defined( $param->{ omit_p_value } )     ? $param->{ omit_p_value }     :     0.10;   # P > 0.1 to mode for omission from he calc
    my $p_val_max_len    = defined( $param->{ p_val_max_len } )    ? $param->{ p_val_max_len }    : undef;      # no limit on length in chi-square
    my $reach_p_value    = defined( $param->{ reach_p_value } )    ? $param->{ reach_p_value }    :     0.05;   # allow direction to drift

    if ( $bias_stats ) { $all_he_modes = 0; $max_iter = 0; $omit_p_value = 1 }

    $log = \*STDERR if ( $log || $log_label ) && ref( $log ) ne 'GLOB';
    $log_label = $genome_title || 'Unidentified DNA' if $log && ! $log_label;

    #-------------------------------------------------------------------------------
    #  Index DNA and translate to amino acids:
    #-------------------------------------------------------------------------------

    $dna && ( ref( $dna ) eq 'ARRAY' ) && @$dna
        or print STDERR "native_codon_usage() called with bad 'dna' parameter value.\n"
           and return ();

    my %dna = map { $_->[0] => $_ } @$dna;
    my $n_gene = @$dna;

    #  If $max_keep is fraction of genes, adjust to the number of genes.

    $max_keep = int( $max_keep * $n_gene + 0.5 ) if $max_keep <= 1;

    my $aa;
    my @aa  = map { $aa = gjoseqlib::translate_seq( $_->[2], 1 );
                    $aa =~ s/\*$//;
                    [ @$_[0,1], $aa ]
                  }
              @$dna;

    #-------------------------------------------------------------------------------
    #  Whole genome mode:
    #-------------------------------------------------------------------------------

    if ( $log )
    {
        print  $log "$log_label\n";
        printf $log "%7d total genes.\n\n", $n_gene;
        print  $log "Finding overall modal codon usage.\n";
    }

    #  Build: [ id, def, counts ]
    my @labeled_counts = map { [ @$_[0,1], gjocodonlib::seq_codon_count_package( $_->[2] ) ] } @$dna;
    #  Just the counts
    my @gen_counts     = map { $_->[2] }  @labeled_counts;
    #  Counts indexed by id
    my %gen_counts     = map { @$_[0,2] } @labeled_counts;
    #  Counts with full description
    my @gen_cnt_ids    = map { [ $_->[2], join( ' ', grep { $_ } @$_[0,1] ) ] } @labeled_counts;

    my $save_cnt_file = $counts_file ? 1 : 0;
    if ( $counts_file )
    {
        $counts_file .= '.counts' if $counts_file !~ /\./;   #  Add extension if necessary
    }
    else
    {
        $counts_file = sprintf( "native_codon_usage_tmp_%09d.counts", int(1e9*rand()) );
    }

    open( CNT, ">$counts_file" )
        or print STDERR "native_codon_usage could not open '$counts_file' for writing.\n"
           and exit;
    if ( $save_cnt_file ) { foreach ( @gen_cnt_ids ) { report_counts( \*CNT, @$_ ) } }
    else                  { foreach ( @gen_counts  ) { report_counts( \*CNT,  $_ ) } }
    close CNT;

    my $mode_opts = { count_file => $counts_file,
                      exponent   => $mode_exponent
                    };
    my $gen_mode = gjocodonlib::modal_codon_usage( \@gen_counts, $mode_opts );

    #-------------------------------------------------------------------------------
    #  Genes matching genome average:
    #-------------------------------------------------------------------------------

    my ( $aver_freq, $n_match_average );
    if ( $average || $log )
    {
        $aver_freq = gjocodonlib::count_to_freq( gjocodonlib::sum_counts( \@gen_counts ), 1 );
        $n_match_average = grep { $_->[1] >= $match_p_value }
                           cnts_and_p_values( $aver_freq, \@gen_counts, $p_val_max_len );
        
    }

    #-------------------------------------------------------------------------------
    #  Genes matching genome mode:
    #-------------------------------------------------------------------------------
    #  Genes matching genome mode (this might have different p-value than omission)

    my @mode_p_vals = cnts_and_p_values( $gen_mode, \@gen_counts, $p_val_max_len );
    my %match_mode = map { $_->[1] >= $match_p_value ? ( $_->[0] => 1 ) : () } @mode_p_vals;
    my $n_match = keys %match_mode;

    #  Record genes to be omitted from high expression calculation due to
    #  similarity to genome mode

    my %omit = map { $_->[1] >= $omit_p_value ? ( $_->[0] => 1 ) : () } @mode_p_vals;
    my $n_omit = keys %omit;

    if ( $log )
    {
        printf $log "%7d genes match (%.1f%%) the genome average (P >= $match_p_value).\n",
                     $n_match_average, 100*$n_match_average/$n_gene;
        printf $log "%7d genes match (%.1f%%) the overall mode (P >= $match_p_value).\n",
                     $n_match, 100*$n_match/$n_gene;
        printf $log "%7d genes marked for omission due to matching the mode (P >= $omit_p_value).\n",
                     $n_omit;
        print  $log "\n";
    }

    #-------------------------------------------------------------------------------
    #  Find homologs of high expression genes in reference genome:
    #-------------------------------------------------------------------------------

    if ( $log ) { print  $log "Initial high expression representatives.\n" }

    my @he_id;
    if ( $initial_he_ids && ref( $initial_he_ids ) eq 'ARRAY' )
    {
        my %seen;
        @he_id = grep { $dna{ $_ } && ! $seen{ $_ }++ } @$initial_he_ids;
        if ( $log )
        {
            printf $log "%7d high expression genes defined by calling program.\n",
                         scalar @$initial_he_ids;
            printf $log "%7d have ids in the DNA data.\n",
                         scalar @he_id;
        }
    }

    if ( @he_id < $min_he_genes )
    {
        $he_ref_module =~ s/\.pm$//;    # remove .pm, if present
        eval( "use $he_ref_module qw( \&ref_aa \&ref_he )" );
        my @ref_aa   = &high_expression_ref::ref_aa();
        my @ref_he   = &high_expression_ref::ref_he();

        my %ref_he   = map { $_ => 1 } @ref_he;
        my $n_ref_he = @ref_he;

        my $bbh_opt = { blast_opts    => $bbh_blast_opts,
                        max_e_value   => $bbh_e_value,
                        min_coverage  => $bbh_coverage,
                        min_identity  => $bbh_identity,
                        min_positives => $bbh_positives,
                        program       => 'blastp',
                        subset        => \@ref_he
                      };

        my ( $bbh ) = bidir_best_hits::bbh( \@ref_aa, \@aa, $bbh_opt );

        @he_id = map { $_->[1] } grep { $ref_he{ $_->[0] } } @$bbh;

        if ( $log )
        {
            printf $log "%7d high expression genes defined in reference genome.\n",
                         $n_ref_he;
            printf $log "%7d bidirectional best hits to high expression genes identified.\n",
                         scalar @he_id;
        }

        $param->{ initial_he_ids } = [ @he_id ];  #  Potentially dicy, but ....
    }

    my $n_orig_he = @he_id;

    if ( $log ) { print  $log "\n" }

    #-------------------------------------------------------------------------------
    #  High expression gene mode:
    #-------------------------------------------------------------------------------

    my @he_counts        = map { $gen_counts{ $_ } } @he_id;
    my %orig_he          = map { $_ => 1 } @he_counts;
    my %orig_nonmodal_he = map { $_ => 1 } grep { ! $match_mode{ $_ } } @he_counts;

    my $pipe = open_evaluation_pipe( $counts_file, $p_val_max_len );
    my $he_mode_opts = {};
    my @he_modes = ();
    my $max_orig_he_cand_match = -1;
    my $rollback = 0;

    my $all    = 0;
    my $done   = 0;
    my $n_iter = 0;
    while ( ! $done )
    {
        #  Remove genes matching mode

        my $n_cand = @he_counts;
        if ( $log )
        {
            printf $log "High expression codon frequencies, estimate $n_iter:\n"
                      . "%7d candidate high expression representatives.\n",
                         $n_cand;
        }

        if ( $n_cand < $min_he_genes )
        {
            if ( $log )
            {
                print $log "\n"
                         . "Too few high expression candidates to compute mode (< $min_he_genes).\n"
                         . "Terminating.\n\n";
            }
            last;
        }

        my @he_nonmodal = grep { ! $omit{ $_ } } @he_counts;
        my $n_survive = @he_nonmodal;
        my $min_needed = max( $min_he_genes, $min_nonmodal * $n_cand );
        if ( $n_survive < $min_needed ) { $all = 1; $done = 1 }
        else                            { @he_counts = @he_nonmodal }

        if ( $log )
        {
            printf $log "%7d candidate high expression genes are distinct from mode (P < $omit_p_value).\n"
                      . "        %.1f%% of high expression candidates are distinct from mode at this P-value.\n",
                         $n_survive, 100*$n_survive/$n_cand;

            print  $log "    *** Computing mode from all candidates and terminating.\n" if $done;
            print  $log "\n";
        }

        #  Modal usage of the candidate high expression gene set:

        my $he_mode = modal_codon_usage( \@he_counts, $he_mode_opts );
        push @he_modes, [ $he_mode, "High expression mode $n_iter" . ( $all ? " (including those matching the mode)" : '' ) ];

        # Match of all genes to the new high expression estimate:

        my @match_he = map { $_->[1] >= $match_p_value ? $_->[0] : () }
                       cnts_and_p_values( $he_mode, \@gen_counts, $p_val_max_len );
        my $n_he = @match_he;

        my $n_both = grep { $match_mode{ $_ } } @match_he;

        #  Match of all genes to the mode -> high expression axis.
        #  Each projection is [ $cnt, $x, $p ]

        my @projections = cnts_x_and_p( $pipe, $gen_mode, $he_mode, \@gen_counts, $p_val_max_len );

        #  Match to high expression axis with x >= 0 (if x < 0, then test is
        #  matching the mode).

        my @match_axis = grep { $match_mode{ $_->[0] } || ( $_->[1] >= 0 && $_->[2] >= $match_p_value ) }
                         @projections;
        my $n_match_axis = @match_axis;

        my $excess_overlap = ( ( $n_both > $max_mode_overlap * $n_he ) ? 0 : 0 );  ### This is inactive for now

        #  Count the original high expression candidates that survive.  We use this
        #  to detect systematic drift away from high expression genes.  If detected,
        #  we rollback the high expression mode to that based on the previous set.

        my $n_orig_he_match          = grep { $orig_he{ $_ } }          map { $_->[0] } @match_axis;
        my $n_orig_he_nonmodal_match = grep { $orig_nonmodal_he{ $_ } } map { $_->[0] } @match_axis;
        $max_orig_he_cand_match = $n_orig_he_nonmodal_match if $n_orig_he_nonmodal_match > $max_orig_he_cand_match;

        my $he_match_loss = ( $n_orig_he_nonmodal_match < ( 1-$max_he_decline ) * $max_orig_he_cand_match );

        $excess_overlap += 2 if ( $n_orig_he_nonmodal_match < ( 1-$max_mode_overlap ) * $n_orig_he_match );

        if ( $log )
        {
            my $prefix = ( $excess_overlap & 1 ) ? '*** Only' : '   ';
            printf $log "%7d genes match high expression estimate $n_iter (P >= $match_p_value).\n"
                      . "%7d genes match both the mode and high expression estimate.\n"
                      . "    $prefix %.1f%% of genes matching the high expression mode are distinct from\n"
                      . "           the overall mode.\n\n",
                         $n_he, $n_both, 100*(1-$n_both/$n_he);

            printf $log "%7d genes (%.1f%%) match the mode - high expression axis (P >= $match_p_value).\n"
                      . "%7d of the original %d high expression candidates match the axis.\n"
                      . "%7d of the original high expression candidates match the axis but not mode.\n",
                         $n_match_axis, 100*$n_match_axis/$n_gene,
                         $n_orig_he_match, $n_orig_he,
                         $n_orig_he_nonmodal_match;

            if ( $excess_overlap & 2 )
            {
                printf $log "    *** Only %.1f%% of the original high expression candidates that match the\n"
                          . "           axis are distinct from the overall mode.\n",
                             100 * $n_orig_he_nonmodal_match/$n_orig_he_match;
            }
            if ( $he_match_loss )
            {
                printf $log "    *** Matches of original high expression candidates has dropped more than %.1f%%.\n",
                             100 * $max_he_decline;
            }
            if ( ( $he_match_loss || $excess_overlap ) && ! $done )
            {
                print  $log "    *** Using first high expression frequencies estimate.\n";
            }
            print  $log "\n";
        }

        $rollback = $excess_overlap || $he_match_loss;
        last if ( $n_iter >= $max_iter ) || $done || $rollback;
        $n_iter++;

        # Limit by p-value (this is generally more relaxed than the match p-value):
        @projections = grep { $_->[2] >= $reach_p_value }  @projections;
        my $pass_p_val = @projections;

        # Limit by x-value (D = 0.5)
        @projections = grep { $_->[1] >= $min_x_value } @projections;
        my $pass_x_val = @projections;

        my $rank = 0;
        @projections = grep { ++$rank <= $max_keep }      # Limit number kept
                       sort { $b->[1] <=> $a->[1] }       # Highest x-value to lowest
                       @projections;
        my $pass_keep = @projections;

        @he_counts = map  { $_->[0] } @projections;

        if ( $log )
        {
            printf $log "New high expression representatives:\n"
                      . "%7d genes pass P-value test (P >= $reach_p_value).\n"
                      . "%7d genes pass x-value test (x >= $min_x_value).\n"
                      . "%7d genes pass max keep test (n <= $max_keep).\n\n",
                         $pass_p_val, $pass_x_val, $pass_keep;
        }
    }

    if ( @he_modes > 1 )
    {
        @he_modes = ( $he_modes[ 0] ) if   $rollback;
        @he_modes = ( $he_modes[-1] ) if ! $all_he_modes;
    }

    close_pipe2( $pipe ) if   $pipe;
    unlink $counts_file  if ! $save_cnt_file;

    if ( $bias_stats && @he_modes )
    {
        #  Find distance between mode and he mode:

        my $dist = gjocodonlib::codon_freq_distance( $gen_mode, $he_modes[0]->[0] );

        #  Do same for 10 random subsets of initial gene set.

        my $set_size = @{ $param->{ initial_he_ids } };
        my @dists;

        for ( my $i = 0; $i < 10; $i++ )
        {
            my @cnts2 = ( sort { rand() <=> 0.5 } @gen_counts )[ 0 .. ($set_size-1) ];
            my $freq2 = gjocodonlib::modal_codon_usage( \@cnts2, {} );
            push @dists, gjocodonlib::codon_freq_distance( $gen_mode, $freq2 );
        }

        my ( $mean, $stddev ) = gjostat::mean_stddev( @dists );
        my $Z = ( $dist - $mean ) / $stddev;
        my $stats = sprintf "distance = %5.3f; random = %5.3f +/- %5.3f; Z = %.1f", $dist, $mean, $stddev, $Z;
        $he_modes[0]->[1] = "High expression mode: $stats";
    }

    return map { $_->[1] =~ s/^/$genome_title -- / if $genome_title; $_ }
           ( ( $average ? [ $aver_freq, "Genome average" ] : () ),
             [ $gen_mode, 'Genome mode' ],
             @he_modes
           );
}


#===============================================================================
#  Just subroutines below:
#===============================================================================
#
#   $pipe = open_evaluation_pipe( $counts_file, $l_max );
#
#   $pipe = [ $pid, $rw, $wr ]
#
#-------------------------------------------------------------------------------
sub open_evaluation_pipe
{
    my ( $counts_file, $l_max ) = @_;
    return () if ! &gjocodonlib::version( 'project_codon_usage_onto_axis'  );

    my ( $rd, $wr, $eval_cmd );
    $eval_cmd  = "project_codon_usage_onto_axis";
    $eval_cmd .= " -l $l_max" if $l_max;
    $eval_cmd .= " $counts_file";
    my $pid = open2( $rd, $wr, $eval_cmd );
    { my $old = select $wr; $| = 1; select $old; }  #  Autoflush the write pipe

    [ $pid, $rd, $wr ];
}

#-------------------------------------------------------------------------------
#  close_pipe2( $pipe );
#-------------------------------------------------------------------------------
sub close_pipe2
{
    my ( $pipe ) = @_;
    if ( $pipe && ( ref( $pipe ) eq 'ARRAY' ) && @$pipe )
    {
        my ( $pid, $rd, $wr ) = @$pipe;
        close $wr if $wr;
        close $rd if $rd;
        waitpid( $pid, 0 ) if $pid;
    }
}

#-------------------------------------------------------------------------------
#
#  @cnt_pval_pairs = cnts_and_p_values( $freq, \@cnts, $l_max );
#  @cnt_pval_pairs = cnts_and_p_values( $freq, \@cnts );
#
#-------------------------------------------------------------------------------
sub cnts_and_p_values
{
    my ( $freq, $cnts, $l_max ) = @_;
    $freq && ref( $freq ) eq 'ARRAY' or return ();
    $cnts && ref( $cnts ) eq 'ARRAY' or return ();
    my ( $chisqr, $df, $n );

    map { ( $chisqr, $df, $n ) = gjocodonlib::count_vs_freq_chi_sqr( $_, $freq );
          $chisqr *= $l_max / $n if $l_max && ( $n > $l_max );
          [ $_, $df ? gjostat::chisqr_prob( $chisqr, $df ) : 1 ]
        } @$cnts;
}

#-------------------------------------------------------------------------------
#
#  ( [ $cnt, $x, $p ], ... ) = cnts_x_and_p( $pipe, $freq1, $freq2, $cnts, $l_max );
#  ( [ $cnt, $x, $p ], ... ) = cnts_x_and_p( $pipe, $freq1, $freq2, $cnts );
#
#-------------------------------------------------------------------------------
sub cnts_x_and_p
{
    my ( $pipe, $freq1, $freq2, $cnts, $l_max ) = @_;
    $freq1 && ref( $freq1 ) eq 'ARRAY' or return ();
    $freq2 && ref( $freq2 ) eq 'ARRAY' or return ();
    $cnts  && ref( $cnts )  eq 'ARRAY' or return ();

    if ( $pipe )
    {
        my ( undef, $rd, $wr ) = @$pipe;
        gjocodonlib::report_frequencies( $wr, $freq1 );
        gjocodonlib::report_frequencies( $wr, $freq2 );
        my $x_p;
        map { chomp( $x_p = <$rd> ); [ $_, split( /\t/, $x_p ) ] } @$cnts;
    }
    else
    {
        my ( $proj, $x, $chisqr, $df, $len );
        map { ( $proj ) = gjocodonlib::project_on_freq_vector_by_chi_sqr_2( $freq1, $freq2, $_ );
              ( $x, $chisqr, $df, $len ) =  @$proj;
              $chisqr *= $l_max / $len if $l_max && ( $len > $l_max );
              [ $_, $x, ( $df ? gjostat::chisqr_prob( $chisqr, $df ) : 1 ) ]
            } @$cnts;
    }
}


#-------------------------------------------------------------------------------
#  $max = max( $n1, $n2 );
#-------------------------------------------------------------------------------
sub max { $_[0] >= $_[1] ? $_[0] : $_[1] }


1;
