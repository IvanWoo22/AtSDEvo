#!/usr/bin/env python
import sys
import argparse
import pandas as pd


def main(input_dats, output_bed):
    header = (
        '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif '
        'Sequence').split()

    trf = []
    for datf in input_dats:
        chrom = None
        sys.stderr.write("\r" + datf)
        with open(datf, 'r') as dat:
            for line in dat:
                splitline = line.split()
                if line.startswith("Sequence:"):
                    chrom = int(line.split()[1].strip())
                elif line.startswith("@"):
                    chrom = splitline[0][1:].strip()
                else:
                    try:
                        int(splitline[0])
                    except ValueError:
                        continue
                    trf.append([chrom] + splitline[0: (len(header) - 1)])

    trf_df = pd.DataFrame(trf, columns=header)
    trf_df["start"] = trf_df["start"].astype(int)
    trf_df.sort_values(by=["#chr", "start"], inplace=True)
    trf_df.to_csv(output_bed, sep="\t", index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TRF output files to BED format.")
    parser.add_argument("input_dats", nargs="+", help="Input TRF data files")
    parser.add_argument("output_bed", help="Output BED file")
    args = parser.parse_args()

    main(args.input_dats, args.output_bed)
