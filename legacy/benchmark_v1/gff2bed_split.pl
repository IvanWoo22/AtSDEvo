#!/usr/bin/env perl
use strict;
use warnings;

my %gene_type;
open my $LIST, "<", $ARGV[0];
while (<$LIST>) {
    chomp;
    my @tmp = split "\t";
    $gene_type{ $tmp[0] } = $tmp[1];
}
open my $GFF, "<", $ARGV[1];
my $BED_SUFFIX = $ARGV[2];
open my $BED_PMT,      ">", $BED_SUFFIX . "/promoter.raw.bed";
open my $BED_GENEBODY, ">", $BED_SUFFIX . "/genebody.raw.bed";
open my $BED_GENE,     ">", $BED_SUFFIX . "/gene.raw.bed";
my %type_check = (
    "novel_transcribed_region"     => 1,
    "long_noncoding_rna"           => 1,
    "antisense_long_noncoding_rna" => 1,
    "protein_coding"               => 1,
    "pseudogene"                   => 1
);

while (<$GFF>) {
    if ( ( !/^#/ ) and ( !/^Chr[MC]/ ) and ( !/^\s*$/ ) ) {
        chomp;
        my @tmp = split "\t";
        if ( ( $tmp[2] eq "gene" ) or ( $tmp[2] eq "pseudogene" ) ) {
            if ( $tmp[8] =~ /ID=([^;]+)(?:;|$)/ ) {
                my $cid = $1;
                if ( exists( $gene_type{$cid} ) ) {
                    my $type = $gene_type{$cid};
                    if ( exists( $type_check{$type} ) ) {
                        $tmp[0] =~ s/Chr//;
                        my $chr = $tmp[0];
                        print $BED_GENEBODY
                          "$chr\t$tmp[3]\t$tmp[4]\t$cid|$type\n";
                        my ( $start, $end );
                        if ( $tmp[6] eq "+" ) {
                            $start = $tmp[3] - 2000;
                            $start = 1 if $start < 1;
                            $end   = $tmp[3] - 1;
                            print $BED_GENE
                              "$chr\t$start\t$tmp[4]\t$cid|$type\n";
                        }
                        else {
                            $start = $tmp[4] + 1;
                            $end   = $tmp[4] + 2000;
                            print $BED_GENE "$chr\t$tmp[3]\t$end\t$cid|$type\n";
                        }
                        print $BED_PMT "$chr\t$start\t$end\t$cid|$type\n";
                    }
                }
            }
        }
    }
}
close($GFF);
close($BED_GENE);
close($BED_PMT);
close($BED_GENEBODY);

__END__
