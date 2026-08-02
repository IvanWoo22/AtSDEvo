#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my (
    $abbr,    $ptag,    $gff_file,  $pep_file,
    $out_gff, $out_pep, %name_info, $out_list
);
Getopt::Long::GetOptions(
    "gff=s"      => \$gff_file,
    "pep=s"      => \$pep_file,
    "abbr=s"     => \$abbr,
    "pep_tag=s"  => \$ptag,
    "out_gff=s"  => \$out_gff,
    "out_pep=s"  => \$out_pep,
    "out_list=s" => \$out_list
  )
  or die
"Usage: $0 --gff GFF_FILE --pep PEP_FILE --abbr ABBR --pep_tag PEP_TAG --out_gff OUT_GFF --out_pep OUT_PEP --out_list OUT_LIST\n";

open my $GFF, "<", $gff_file or die "Cannot open GFF file: $!";
open my $OGF, ">", $out_gff  or die "Cannot create output GFF file: $!";
while (<$GFF>) {
    chomp;
    my ($name) = /$ptag([^;]+)(?:;|\t)/;
    my @tmp    = split( /\t/, $_ );
    my $contig = $tmp[0];
    my ( $start, $end ) = ( $tmp[1], $tmp[2] );
    ( $start, $end ) = ( $end, $start ) if $tmp[4] eq "-";
    $name_info{$name} = "$abbr\_$contig\t$name\t$start\t$end\n";
}
close $GFF;

open my $PEP,   "<", $pep_file or die "Cannot open PEP file: $!";
open my $OPEP,  ">", $out_pep  or die "Cannot create output PEP file: $!";
open my $OLIST, ">", $out_list or die "Cannot create output LIST file: $!";
my ( %pep_hash, $current_name, $add );
while (<$PEP>) {
    chomp;
    if (/^>(.+?)(?:\.p)?\s/) {
        $current_name = $1;
        if ( exists $pep_hash{$current_name} ) {
            $add = 0;
        }
        else {
            $add = 1;
            $pep_hash{$current_name} = "";
            print $OLIST "$current_name\n";
        }
    }
    else {
        if ( $add != 0 ) {
            $pep_hash{$current_name} .= $_;
        }
    }
}
foreach my $name ( keys %name_info ) {
    if ( exists $pep_hash{$name} ) {
        print $OGF "$name_info{$name}";
        print $OPEP ">$name\n$pep_hash{$name}\n";
    }
}
close $OGF;
close $PEP;
close $OPEP;
close $OLIST;
