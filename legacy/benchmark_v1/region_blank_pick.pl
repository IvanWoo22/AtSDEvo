#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'fa|f=s'     => \my $in_fa,
    'region|r=s' => \my $in_region,
) or Getopt::Long::HelpMessage(1);

sub FAS2ARR {
    my ( $QUERY, $START, $END, $TARGET ) = @_;
    my @ARRAY1 = split( "", $QUERY );
    my @ARRAY_ABS;
    foreach my $KEY ( 0 .. $#ARRAY1 ) {
        if ( $ARRAY1[$KEY] ne "-" ) {
            push( @ARRAY_ABS, $KEY );
        }
    }
    my @ARRAY2 = split( //, $TARGET );
    return (
        join( "", @ARRAY1[ @ARRAY_ABS[ $START - 1 .. $END - 1 ] ] ),
        join( "", @ARRAY2[ @ARRAY_ABS[ $START - 1 .. $END - 1 ] ] )
    );
}

sub CALCULATE_DASH_RATIO {
    my ($SEQ)        = @_;
    my $TOTAL_LENGTH = length($SEQ);
    my $DASH_COUNT   = ( $SEQ =~ tr/-// );
    my $RATIO        = $TOTAL_LENGTH > 0 ? $DASH_COUNT / $TOTAL_LENGTH : 0;
    return $RATIO;
}

my %fasta;
my %pair;
my @pair_name;
my $title_name;

open my $FA, "<", $in_fa;
while (<$FA>) {
    $_ =~ s/\r?\n//;
    if (/^>(\S+)/) {
        $title_name = $1;
        push( @pair_name, $title_name );
    }
    elsif (/\S+/) {
        $fasta{$title_name} .= $_;
    }
    elsif ( $_ eq '' ) {
        $pair{ $pair_name[0] } = $pair_name[1];
        $pair{ $pair_name[1] } = $pair_name[0];
        @pair_name             = ();
    }
}
close($FA);

open my $RG, "<", $in_region;
while (<$RG>) {
    $_ =~ s/\r?\n//;
    my @tmp = split( "\t", $_ );
    my ( $query_seq, $target_seq ) = FAS2ARR( $fasta{ $tmp[0] },
        $tmp[1], $tmp[2], $fasta{ $pair{ $tmp[0] } } );
    my $ratio1 = CALCULATE_DASH_RATIO($query_seq);
    my $ratio2 = CALCULATE_DASH_RATIO($target_seq);
    print
"$tmp[7]\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$ratio1\t$ratio2\n";
}
close($RG);

__END__
