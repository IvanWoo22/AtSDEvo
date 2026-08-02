#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

open IN, "<", $ARGV[0];
open LK, "<", $ARGV[1];
open OT, ">", $ARGV[2];

my @sit_chr;
my %sit_meth;
foreach ( 1 .. 5 ) {
    $sit_chr[$_] = AlignDB::IntSpan->new;
}

while (<IN>) {
    chomp;
    my @tmp = split( /\t/, $_ );
    $sit_chr[ $tmp[0] ]->AlignDB::IntSpan::add( $tmp[1] );
    my $tmp = $tmp[0] . "S" . $tmp[1];
    $sit_meth{$tmp} = $tmp[3];
}
close IN;

my $line_number = $ARGV[3];
while (<LK>) {
    chomp;
    my @lin = split( "\t", $_ );
    my $set = AlignDB::IntSpan->new;
    $set->add( $lin[ $line_number - 1 ] );
    my $sum      = $sit_chr[ $lin[0] ]->AlignDB::IntSpan::overlap($set);
    my $arg_meth = "NULL";

    if ( $sum > 0 ) {
        my @intsec =
          $sit_chr[ $lin[0] ]->AlignDB::IntSpan::intersect($set)->as_array;
        my $meth = 0;
        foreach ( 0 .. $#intsec ) {
            my $tmp = $lin[0] . "S" . $intsec[$_];
            $meth += $sit_meth{$tmp};
        }
        $arg_meth = $meth / ( $#intsec + 1 );
    }
    print OT "$_\t$sum\t$arg_meth\n";
}
close LK;
close OT;

__END__
