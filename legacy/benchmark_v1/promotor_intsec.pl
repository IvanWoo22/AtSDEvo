#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use AlignDB::IntSpan;

# Variables for command-line options
my $ref_bed;
my $in_bed;
my $output;
my $diff = 0;

# Get options from command line
GetOptions(
    'ref_bed=s' => \$ref_bed,
    'in_bed=s'  => \$in_bed,
    'output=s'  => \$output,
    'diff'      => \$diff,
) or die "Error in command line arguments";

# Check if all required arguments are provided
die "Usage: $0 --ref_bed <ref.bed> --in_bed <in.bed> --output <output.bed>\n"
  unless $ref_bed
  and $in_bed
  and $output;

# Create IntSpan objects for chromosomes
my @chr_set;
foreach ( 1 .. 5 ) {
    $chr_set[$_] = AlignDB::IntSpan->new;
}

# Read the first input file
open my $REF_BED_FH, '<', $ref_bed or die "Cannot open $ref_bed: $!";
while (<$REF_BED_FH>) {
    chomp;
    my @tmp   = split( "\t", $_ );
    my $chr   = $tmp[0];
    my $start = $tmp[1];
    my $end   = $tmp[2];
    $chr_set[$chr]->AlignDB::IntSpan::add_pair( $start, $end );
}
close $REF_BED_FH;

# Read the second input file
open my $IN_BED_FH, '<', $in_bed or die "Cannot open $in_bed: $!";
open my $OUT_FH,    '>', $output or die "Cannot open $output: $!";
while (<$IN_BED_FH>) {
    chomp;
    my @tmp = split( "\t", $_ );
    my $set = AlignDB::IntSpan->new;
    $set->add_pair( $tmp[1], $tmp[2] );
    my $int_set;
    if ( $diff != 0 ) {
        $int_set = $set->diff( $chr_set[ $tmp[0] ] );
    }
    else {
        $int_set = $set->intersect( $chr_set[ $tmp[0] ] );
    }
    if ( $int_set->cardinality >= 5 ) {
        my $rl = $int_set->as_string;
        print $OUT_FH "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$rl\n";
    }
}
close $IN_BED_FH;
close $OUT_FH;

__END__
