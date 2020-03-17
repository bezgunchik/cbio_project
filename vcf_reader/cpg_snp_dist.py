import sys
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import math

def create_matrix(input_path):
    file = open(input_path, "r")
    cpg_dict = {}
    snp_dict = {}

    # line = file.readline()
    for line in file:
        split_line = line.split(',')
        cpg_dict[split_line[0]] = split_line[1]
        snp_dict[split_line[2]] = split_line[3].split("\n")[0]
    file.close()
    # matrix = []
    # row = []
    # row.append("")
    # for cpg in cpg_dict:
    #     row.append(cpg)
    # matrix.append(row)
    #
    # for snp in snp_dict:
    #     row = []
    #     row.append(snp)
    #     for cpg in cpg_dict:
    #         row.append(str(abs(int(snp_dict.get(snp)) - int(cpg_dict.get(cpg)))))
    #     matrix.append(row)
    # print(matrix)

    snp_names = []
    cpg_names = []

    for cpg in cpg_dict:
        cpg_names.append(cpg)


    matrix = []
    row = []


    for snp in snp_dict:
        snp_names.append(snp)
        row = []
        for cpg in cpg_dict:
            row.append(math.log10(abs(int(snp_dict.get(snp)) - int(cpg_dict.get(cpg))) + 1))
        matrix.append(row)


    # table = pd.DataFrame(matrix, index=snp_names, columns=cpg_names)
    table = pd.DataFrame(matrix, index=snp_names, columns=cpg_names)


    print(table)
    sb.set(font_scale=1.5)

    sb.heatmap(table, xticklabels=False, yticklabels=False, cbar_kws={'label': 'log10'})
    # plt.rcParams['font.size'] = '50'
    plt.xlabel("CpG", {'fontname':'Calibri'} )
    plt.ylabel("SNP", {'fontname':'Calibri'} )
    plt.title("Distance between chosen SNPs and CpGs", {'fontname':'Calibri'})
    plt.figure()
    plt.show()
    print(len(snp_names), len(cpg_names))
    # print(cpg_dict)
    # print(snp_dict)
# print(file.readline())
    # while line != None:
    #     line.
    #     cpg_list
    # pd.read_csv()


if __name__ == '__main__':
    input_path = sys.argv[1]
    create_matrix(input_path)