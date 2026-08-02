#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

my $partial = $ARGV[-1];

sub READ_YAML {
    my $fh  = shift;
    my @raw = <$fh>;
    my @out;
    foreach ( 1 .. $#raw ) {
        chomp( $raw[$_] );
        my @tmp = split( ": ", $raw[$_] );
        $out[$_] = AlignDB::IntSpan->new;
        $out[$_]->AlignDB::IntSpan::add_runlist( $tmp[1] );
    }
    return ( \@out );
}

sub COVER_JUDGE {
    my ( $set1, $set2, $chr1, $chr2, $yml ) = @_;
    my $comp;
    if (
        (
            $set1->AlignDB::IntSpan::overlap( ${$yml}[$chr1] ) >=
            $partial * $set1->AlignDB::IntSpan::size
        )
        and ( $set2->AlignDB::IntSpan::overlap( ${$yml}[$chr2] ) >=
            $partial * $set2->AlignDB::IntSpan::size )
      )
    {
        $comp = 3;
    }
    elsif ( $set2->AlignDB::IntSpan::overlap( ${$yml}[$chr2] ) >=
        $partial * $set2->AlignDB::IntSpan::size )
    {
        $comp = 2;
    }
    elsif ( $set1->AlignDB::IntSpan::overlap( ${$yml}[$chr1] ) >=
        $partial * $set1->AlignDB::IntSpan::size )
    {
        $comp = 1;
    }
    else {
        $comp = 0;
    }
    return ($comp);
}

my @rl;
foreach ( 1 .. ( $#ARGV - 1 ) ) {
    open my $YAML, "<", $ARGV[$_];
    push( @rl, READ_YAML($YAML) );
    close($YAML);
}

open my $LINKS, "<", $ARGV[0];
while (<$LINKS>) {
    chomp;
    my @tmp  = split( " ", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    my $out1;
    my $out2;
    if ( $tmp[1] < $tmp[2] ) {
        $out1 = $tmp[0] . " " . $tmp[1] . " " . $tmp[2];
        $set1->add_pair( $tmp[1], $tmp[2] );
    }
    else {
        $out1 = $tmp[0] . " " . $tmp[2] . " " . $tmp[1];
        $set1->add_pair( $tmp[2], $tmp[1] );
    }

    if ( $tmp[4] < $tmp[5] ) {
        $out2 = $tmp[3] . " " . $tmp[4] . " " . $tmp[5];
        $set2->add_pair( $tmp[4], $tmp[5] );
    }
    else {
        $out2 = $tmp[3] . " " . $tmp[5] . " " . $tmp[4];
        $set2->add_pair( $tmp[5], $tmp[4] );
    }

    print "$out1\t$out2";
    foreach (@rl) {
        my $out_comp = COVER_JUDGE( $set1, $set2, $tmp[0], $tmp[3], \@{$_} );
        print "\t$out_comp";
    }
    print "\t$tmp[6]\n";
}
close($LINKS);
__END__
