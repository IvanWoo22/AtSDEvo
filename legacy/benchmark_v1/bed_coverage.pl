#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

while (<STDIN>) {
    chomp;
    my @tmp  = split "\t";
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    $set1->add_range( $tmp[1], $tmp[2] );
    $set2->add_range( $tmp[5], $tmp[6] );
    if ( $set1->overlap($set2) / $set1->size > 0.5 ) {
        print( join( "\t", @tmp ) );
        print "\t" . $set1->overlap($set2) / $set1->size . "\n";
    }
}

__END__