#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub OVERLAP_CLOSE {
    my ( $start1, $end1, $start2, $end2, $threshold ) = @_;

    my $length1 = $end1 - $start1 + 1;
    my $length2 = $end2 - $start2 + 1;

    if ( $end1 >= $start2 && $end2 >= $start1 ) {
        my $overlap_start  = $start1 > $start2 ? $start1 : $start2;
        my $overlap_end    = $end1 < $end2     ? $end1   : $end2;
        my $overlap_length = $overlap_end - $overlap_start + 1;

        my $proportion1 = $overlap_length / $length1;
        my $proportion2 = $overlap_length / $length2;

        return ( $proportion1 > $threshold || $proportion2 > $threshold );
    }
    elsif ( $threshold < 0 ) {
        my $distance =
          ( $start2 > $end1 ) ? $start2 - $end1 - 1 : $start1 - $end2 - 1;
        my $proportion1 = $distance / $length1;
        my $proportion2 = $distance / $length2;

        return ( $proportion1 < -$threshold || $proportion2 < -$threshold );
    }

    return 0;
}

my $input_file;
my $overlap_threshold = 0.2;
my $exclude           = 0;

GetOptions(
    'input=s'   => \$input_file,
    'overlap=f' => \$overlap_threshold,
    'exclude'   => \$exclude,
) or die "Error in command line arguments\n";

if ( !defined $input_file ) {
    die "Usage: $0 --input <input_file> [--overlap <proportion>] [--exclude]\n";
}

open my $FH, '<', $input_file or die "Cannot open file $input_file: $!";

while ( my $line = <$FH> ) {
    chomp($line);
    my @intervals               = split( /\s+/, $line );
    my $has_significant_overlap = 0;

    for ( my $i = 0 ; $i < @intervals ; $i++ ) {
        my ( $chr1,   $pos1 ) = split( /:/, $intervals[$i] );
        my ( $start1, $end1 ) = split( /-/, $pos1 );
        for ( my $j = $i + 1 ; $j < @intervals ; $j++ ) {
            my ( $chr2,   $pos2 ) = split( /:/, $intervals[$j] );
            my ( $start2, $end2 ) = split( /-/, $pos2 );
            if (
                $chr1 eq $chr2
                && OVERLAP_CLOSE(
                    $start1, $end1, $start2, $end2, $overlap_threshold
                )
              )
            {
                $has_significant_overlap = 1;
                last;
            }
        }
        last if $has_significant_overlap != 0;
    }
    if (   ( $exclude && !$has_significant_overlap )
        || ( !$exclude && $has_significant_overlap ) )
    {
        print "$line\n";
    }
}
close($FH);

__END__
