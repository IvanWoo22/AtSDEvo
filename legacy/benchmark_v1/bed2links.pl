#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my @bed_text;
my @dup_bool;
my @other_info;
open my $BED, "<", $ARGV[0];
while (<$BED>) {
    chomp;
    my @tmp = split( "\t", $_ );
    my $dir = "+";
    if ( $tmp[4] != 0 ) {
        if ( ( $tmp[5] ne "F" ) && ( $tmp[4] == 2 ) ) {
            $dir = "-";
        }
        ${ $bed_text[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] =
          "Chr$tmp[0]\t$tmp[1]\t$tmp[2]\t\t1\t$dir";
        ${ $dup_bool[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] = 1;
        ${ $other_info[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] =
          join( "\t", @tmp[ 0, 1, 2, 5, 6 ] );
    }
    else {
        if (   ( $tmp[5] ne "F" )
            && defined( $bed_text[ $tmp[3] - 1 ] )
            && ( @{ $bed_text[ $tmp[3] - 1 ] } ) )
        {
            $dir = "-";
        }
        push(
            @{ $bed_text[ $tmp[3] - 1 ] },
            "Chr$tmp[0]\t$tmp[1]\t$tmp[2]\t\t1\t$dir"
        );
        push( @{ $dup_bool[ $tmp[3] - 1 ] }, 0 );
        push(
            @{ $other_info[ $tmp[3] - 1 ] },
            join( "\t", @tmp[ 0, 1, 2, 5, 6 ] )
        );
    }
}
close($BED);

open my $INFO, ">", $ARGV[2];
foreach my $key ( 0 .. $#bed_text ) {
    if ( defined( $bed_text[$key] ) ) {
        open my $TMP, ">", "$ARGV[3]/$key" . ".bed";
        print $TMP ( join( "\n", @{ $bed_text[$key] } ) );
        print $TMP ("\n");
        print $INFO (
            join(
                "\n",
                map {
                        ( $key + 1 ) . "-"
                      . ( $_ + 1 )
                      . "\t${ $bed_text[$key] }[$_]\t${ $dup_bool[$key] }[$_]\t${$other_info[$key] }[$_]"
                } 0 .. $#{ $bed_text[$key] }
            )
        );
        print $INFO ("\n");
        my $command =
            "bedtools getfasta -s -fi $ARGV[1] -bed $ARGV[3]/$key"
          . ".bed >$ARGV[3]/$key" . ".fa";
        my $status = system($command);
        if ( $status != 0 ) {
            print "Command1 failed with status: $status\n";
        }
        $command = "sed -i 's/([+-])\$//' $ARGV[3]/$key" . ".fa";
        $status  = system($command);
        if ( $status != 0 ) {
            print "Command2 failed with status: $status\n";
        }
        $command = "sed -i 's/:/(+):/' $ARGV[3]/$key" . ".fa";
        $status  = system($command);
        if ( $status != 0 ) {
            print "Command3 failed with status: $status\n";
        }
    }
}
__END__
