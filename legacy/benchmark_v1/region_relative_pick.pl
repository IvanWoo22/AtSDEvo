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

sub FASTA_RELATIVE_REGION {
    my ( $QUERY, $START, $END, $TARGET ) = @_;
    my @ARRAY1 = split( "", $QUERY );
    my @ARRAY_ABS1;
    foreach my $KEY ( 0 .. $#ARRAY1 ) {
        if ( $ARRAY1[$KEY] ne "-" ) {
            push( @ARRAY_ABS1, $KEY );
        }
    }
    my @ARRAY2 = split( "", $TARGET );
    my @ARRAY_ABS2;
    my $KEY = 0;
    foreach ( 0 .. $#ARRAY2 ) {
        push( @ARRAY_ABS2, $KEY );
        if ( $ARRAY2[$_] ne "-" ) {
            $KEY++;
        }
    }
    return (
        $ARRAY_ABS2[ $ARRAY_ABS1[$START] ],
        $ARRAY_ABS2[ $ARRAY_ABS1[ $END - 1 ] ] - 1
    );
}

sub EXTRACT_INFO {
    my ($TEXT) = @_;
    if ( $TEXT =~ /Chr(\d+)\(\+\):(\d+)-\d+/ ) {
        my $CHR   = $1;
        my $START = $2;
        return ( $CHR, $START );
    }
    else {
        return undef;
    }
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
    my ( $target_start, $target_end ) =
      FASTA_RELATIVE_REGION( $fasta{ $tmp[0] },
        $tmp[1], $tmp[2], $fasta{ $pair{ $tmp[0] } } );
    my ( $chr, $start_ori ) = EXTRACT_INFO( $pair{ $tmp[0] } );
    my $abs_start = $start_ori + $target_start;
    my $abs_end   = $start_ori + $target_end;
    print
"$chr\t$abs_start\t$abs_end\t$tmp[7],$tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4],$pair{$tmp[0]},$target_start,$target_end\n";
}
close($RG);

__END__
