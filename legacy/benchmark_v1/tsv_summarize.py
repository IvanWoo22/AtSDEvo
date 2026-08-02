#!/usr/bin/env python
import sys


def summarize_tsv(tsv_file, columns):
    with open(tsv_file, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            count_0 = count_1 = count_2 = count_3 = 0
            for col in columns:
                value = int(data[col - 1])
                if value == 0:
                    count_0 += 1
                elif value == 1:
                    count_1 += 1
                elif value == 2:
                    count_2 += 1
                elif value == 3:
                    count_3 += 1
            if count_3 > 0:
                print(3)
            elif count_2 > 0 and count_1 == 0:
                print(2)
            elif count_1 > 0 and count_2 == 0:
                print(1)
            elif count_1 > 0 and count_2 > 0:
                print(3)
            else:
                print(0)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <TSV_FILE> <COLUMNS>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    columns = [int(col) for col in sys.argv[2].split(",")]

    summarize_tsv(tsv_file, columns)
