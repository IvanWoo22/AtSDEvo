#!/usr/bin/perl
use strict;
use warnings;
use autodie;
my %info;
while (<STDIN>) {
    chomp;
    my @tmp1 = split "\t";
    my ( $gene1, undef ) = split( /\|/, $tmp1[7] );
    my @tmp2 = split( ",", $tmp1[3] );
    my ( $gene2, undef ) = split( /\|/, $tmp2[0] );
    $info{$gene1} = $tmp2[6];
    $info{$gene2} = $tmp2[1];
    my ( $sorted_gene1, $sorted_gene2 ) = sort( $gene1, $gene2 );
    my $paired = $tmp2[5];

    if ( $sorted_gene1 eq $gene1 ) {
        if ( $tmp2[5] == 2 ) {
            $paired = 1;
        }
        elsif ( $tmp2[5] == 1 ) {
            $paired = 2;
        }
        else {
            $paired = 0;
        }
    }
    print(
"$sorted_gene1\t$sorted_gene2\t$info{$sorted_gene1}\t$info{$sorted_gene2}\t$tmp2[4],$paired\n"
    );
}

__END__
