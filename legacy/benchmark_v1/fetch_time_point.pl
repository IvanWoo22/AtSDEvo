#!/usr/bin/env perl
use strict;
use warnings;

my %feature;
open( my $TIME, "<", $ARGV[0] );
while (<$TIME>) {
    chomp;
    my @tmp = split "\t";
    $feature{ $tmp[0] } = $tmp[1];
}
close($TIME);

open( my $LINE, "<", $ARGV[1] );
while (<$LINE>) {
    chomp;
    my @tmp  = split "\t";
    my $time = $tmp[2] . $tmp[3] . $tmp[4] . $tmp[5] . $tmp[6];
    if ( exists( $feature{$time} ) ) {
        if ( $feature{$time} == 0 ) {
            print "$tmp[0]\t$tmp[1]\t0\t$time\n";
        }
        elsif ( $feature{$time} == 1 ) {
            print "$tmp[0]\t$tmp[1]\t1\t$time\n";
        }
        elsif ( $feature{$time} == 2 ) {
            print "$tmp[1]\t$tmp[0]\t1\t$time\n";
        }
    }
}
close($LINE);

__END__
