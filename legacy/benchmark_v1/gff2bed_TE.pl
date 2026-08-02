#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open my $GFF, "<", $ARGV[0];
while (<$GFF>) {
    if ( ( !/^#/ ) and ( !/^Chr[MC]/ ) and ( !/^\s*$/ ) ) {
        chomp;
        my @tmp = split "\t";
        if ( $tmp[2] eq "transposable_element" ) {
            if ( $tmp[8] =~ /ID=([^;]+)(?:;|$)/ ) {
                my $cid  = $1;
                my $type = $tmp[2];
                $tmp[0] =~ s/Chr//;
                my $chr = $tmp[0];
                print "$chr\t$tmp[3]\t$tmp[4]\t$cid|$type\n";
            }
        }
    }
}
close($GFF);

__END__
