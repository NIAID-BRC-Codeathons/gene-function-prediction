#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($min) = @ARGV;
$min ||= 0.0;

my $line;
my %links;
while ($line = <STDIN>) {
    if ($line =~ m/^(\d+\.\d+)/) {
	my $genome = $1;
	print STDERR ($genome, "\n");
	
	# 0: 'genome.genome_id',
	# 1: 'feature.patric_id',
        # 2: 'feature.feature_type',
        # 3: 'feature.sequence_id',
        # 4: 'feature.start',
        # 5: 'feature.end',
        # 6: 'feature.strand',
        # 7: 'feature.product',
        # 8: 'feature.pgfam_id'
	my @out = map { chomp;
			[ split /\t/ ]
	} `p3-echo -t genome.genome_id $genome | p3-get-genome-features --attr patric_id,feature_type,sequence_id,start,end,strand,product,pgfam_id | p3-sort feature.sequence_id feature.start/n feature.end/nr`;
	#die Dumper($out[0]);
	shift @out;   #...throw away header line...
	#die Dumper(\@out);
	
	for (my $i=1; $i < @out; ++$i) {
	    my (undef, $fid_1, $type_1, $contig_1, $left_1, $right_1, $strand_1, $func_1, $pgfam_id_1) = @ { $out[$i-1] };
	    my (undef, $fid_2, $type_2, $contig_2, $left_2, $right_2, $strand_2, $func_2, $pgfam_id_2) = @ { $out[$i] };
	    
	    next unless ($contig_1 eq $contig_2);
	    #next unless ($strand_1 eq $strand_2);
	    next unless ($type_1   eq 'CDS');
	    next unless ($type_2   eq 'CDS');
	    #next unless $pgfam_id_1;
	    #next unless $pgfam_id_2;
	    
	    my $gap = $left_2 - $right_1;
	    next if ($gap > 45);
	    
	    # my ($upfam, $downfam);
	    # if ($strand_1 eq '+') {
	    # 	$upfam   = $pgfam_id_1;
	    # 	$downfam = $pgfam_id_2;
	    # }
	    # else {
	    # 	$upfam   = $pgfam_id_2;
	    # 	$downfam = $pgfam_id_1;
	    # }

	    my @roles_1 = &roles_of_function($func_1);
	    my @roles_2 = &roles_of_function($func_2);

	    
	    foreach my $role_1 (@roles_1) {
		foreach my $role_2 (@roles_2) {
		    my $key = ($role_1 lt $role_2) ? "$role_1\t$role_2" : "$role_2\t$role_1";
		    ++$links{$key};
		}
	    }
	}
    }
}

foreach my $key (sort { $links{$b} <=> $links{$a} } keys %links) {
    print STDOUT ($links{$key}, "\t", $key, "\n");
}
exit(0);


# foreach my $source (sort keys %links) {
#     foreach my $target (sort keys % { $links{$source} } ) {
# 	my $score = $links{$source}->{$target};
# 	if ($score >= $min) {
# 	    print STDOUT (join("\t", ($source, $target, $score)), "\n");
# 	}
#     }
# }

sub roles_of_function {
    # Get the parameters.
    my ($assignment) = @_;
    # Remove any comment.
    my $commentFree = ($assignment =~ /(.+?)\s*[#!]/ ? $1 : $assignment);
    # Split out the roles.
    my @retVal = grep { $_ } split /\s+[\/@]\s+|\s*;\s+/, $commentFree;
    # Return the result.
    return @retVal;
}
