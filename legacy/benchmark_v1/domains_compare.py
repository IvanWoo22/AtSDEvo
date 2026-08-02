#!/usr/bin/env python
import argparse
import ast
from decimal import Decimal, getcontext

getcontext().prec = 4


def read_input(filename):
    protein_data = {}
    with open(filename, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            protein_name = data[0]
            domains = ast.literal_eval(data[1])
            protein_data[protein_name] = domains
    return protein_data


def read_domain_scores(filename):
    domain_scores = {}
    with open(filename, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            if len(data) == 2:
                domain_pair, score = data
                domain_scores[domain_pair] = float(score)
    return domain_scores


def needleman_wunsch(list1, list2, domain_scores):
    len_list1 = len(list1)
    len_list2 = len(list2)

    dp_matrix = [[Decimal("0")] * (len_list2 + 1) for _ in range(len_list1 + 1)]
    gap = Decimal("0")

    for i in range(len_list1 + 1):
        dp_matrix[i][0] = i * gap
    for j in range(len_list2 + 1):
        dp_matrix[0][j] = j * gap
    for i in range(1, len_list1 + 1):
        for j in range(1, len_list2 + 1):
            match_str1 = f"{list1[i - 1]}={list2[j - 1]}"
            match_str2 = f"{list2[j - 1]}={list1[i - 1]}"
            # match_mismatch=0
            if match_str1 in domain_scores:
                match_mismatch_score = Decimal(domain_scores[match_str1])
            elif match_str2 in domain_scores:
                match_mismatch_score = Decimal(domain_scores[match_str2])
            else:
                match_mismatch_score = Decimal("-0.05")
            dp_matrix[i][j] = max(
                dp_matrix[i - 1][j - 1] + match_mismatch_score,
                dp_matrix[i - 1][j] + gap,
                dp_matrix[i][j - 1] + gap
            )
    score = dp_matrix[len_list1][len_list2]

    alignment_list1 = []
    alignment_list2 = []
    i, j = len_list1, len_list2
    while i > 0 or j > 0:
        if i > 0 and dp_matrix[i][j] == dp_matrix[i - 1][j] + gap:
            alignment_list1.insert(0, list1[i - 1])
            alignment_list2.insert(0, "-")
            i -= 1
        elif j > 0 and dp_matrix[i][j] == dp_matrix[i][j - 1] + gap:
            alignment_list1.insert(0, "-")
            alignment_list2.insert(0, list2[j - 1])
            j -= 1
        else:
            alignment_list1.insert(0, list1[i - 1])
            alignment_list2.insert(0, list2[j - 1])
            i -= 1
            j -= 1
    if len(alignment_list1) != 0:
        count = 0
        for i in range(len(alignment_list1)):
            if alignment_list1[i] != "-" and alignment_list2[i] != "-":
                count += 1
        score = Decimal(score) / Decimal(len(alignment_list1))
        # print("match"+ str(non_gap_count))
        # print("list1"+ str(len(alignment_list1)))
        # print("list2"+ str(len(alignment_list2)))
        # print("match"+ str(count))
        score = score * ((Decimal(count) / Decimal(len(alignment_list1))) ** 2)
        # Jaccard相似性
    else:
        score = "NULL"
    return alignment_list1, alignment_list2, score
    # return alignment_list1,score


def main():
    parser = argparse.ArgumentParser(
        description="python protein_domain_blast.py -i input.txt   -I  input_2.txt   -s score.txt -o output.txt ",
        epilog="""输入文件为tsv文件 ： 蛋白名称 domain1 domain2  打分规则文件为tsv文件 ： domain1=domain2 分数（0~1）"""
    )

    parser.add_argument("-i", "--input", help="输入文件1的路径")

    parser.add_argument("-I", "--input_2", help="输入文件2的路径")

    parser.add_argument("-s", "--score", help="打分规则文件")

    parser.add_argument("-o", "--output", help="结果文件")

    # 解析命令行参数
    args = parser.parse_args()
    protein_file1 = args.input
    protein_file2 = args.input_2
    score_file = args.score
    output_file = args.output

    protein_data1 = read_input(protein_file1)
    protein_data2 = read_input(protein_file2)
    domain_scores = read_domain_scores(score_file)
    with open(output_file, "a") as file:
        for protein1_name, protein1 in protein_data1.items():
            for protein2_name, protein2 in protein_data2.items():
                if protein1_name != protein2_name:
                    result_alignment_list1, result_alignment_list2, score = needleman_wunsch(protein1, protein2,
                                                                                             domain_scores)
                    file.write("{}\t{}\t{}\n".format(protein1_name, protein2_name, score))
                    # print(f"{protein1_name}\t{protein2_name}\t{score}")
                    # alignment_str1 = ' '.join(result_alignment_list1)
                    # alignment_str2 = ' '.join(result_alignment_list2)
                    # print(f"\t{alignment_str1}")
                    # print(f"\t{alignment_str2}\n")


if __name__ == "__main__":
    main()
