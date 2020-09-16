import sys
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np
import os

def create_diagram(input_path):
    file = open(input_path)
    snp_dict = {}

    for line in file:
        split_line = line.split(",")
        snp_dict[split_line[0]] = split_line[2].split("\n")[0]
    # counter += 1
    # if counter == 10:
    #     break
    file.close()

    top_snp = {}
    top_snp['AT'] = 0
    top_snp['AC'] = 0
    top_snp['AG'] = 0
    top_snp['TA'] = 0
    top_snp['TC'] = 0
    top_snp['TG'] = 0
    top_snp['CA'] = 0
    top_snp['CT'] = 0
    top_snp['CG'] = 0
    top_snp['GA'] = 0
    top_snp['GT'] = 0
    top_snp['GC'] = 0

    for snp in snp_dict:
        if snp_dict.get(snp) == 'AT':
            top_snp['AT'] += 1
        if snp_dict.get(snp) == 'AC':
            top_snp['AC'] += 1
        if snp_dict.get(snp) == 'AG':
            top_snp['AG'] += 1
        if snp_dict.get(snp) == 'TA':
            top_snp['TA'] += 1
        if snp_dict.get(snp) == 'TC':
            top_snp['TC'] += 1
        if snp_dict.get(snp) == 'TG':
            top_snp['TG'] += 1
        if snp_dict.get(snp) == 'CA':
            top_snp['CA'] += 1
        if snp_dict.get(snp) == 'CT':
            top_snp['CT'] += 1
        if snp_dict.get(snp) == 'CG':
            top_snp['CG'] += 1
        if snp_dict.get(snp) == 'GA':
            top_snp['GA'] += 1
        if snp_dict.get(snp) == 'GT':
            top_snp['GT'] += 1
        if snp_dict.get(snp) == 'GC':
            top_snp['GC'] += 1

    labels = top_snp.keys()
    sizes = top_snp.values()


    # fig = plt.figure()
    # ax = fig.add_axes([0,0,1,1])
    # ax.bar(labels, sizes)

    # counter = 0
    # for line in file:
    #     split_line = line.split(",")
    #     snp_dict[split_line[0]] = float(split_line[1].split("\n")[0])
    #     # counter += 1
    #     # if counter == 10:
    #     #     break
    # file.close()


    # # height =
    # x = snp_dict.keys()
    # # x_pos = [i for i, _ in enumerate(x)]
    # # x = np.array(snp_dict.keys())
    # # y = np.array(snp_dict.values())
    y_pos = np.arange(len(labels))
    plt.bar(labels,sizes, align='center', alpha=0.5)
    plt.ylabel('Number of SNPs')
    plt.xlabel('SNP type')


    plt.title('SNPs near CpGs (TOP 523)')
    # plt.title(os.path.basename(input_path)[11:-4])



    # # plt.xticks(y_pos, snp_dict.keys())
    # # table = pd.DataFrame(snp_dict.values(), columns=snp_dict.keys())
    #
    # # print(table)
    plt.show()

def create_pie(input_path):
    file = open(input_path)
    snp_dict = {}
    counter = 0
    for line in file:
        split_line = line.split(",")
        snp_dict[split_line[0]] = float(split_line[1].split("\n")[0])
        # counter += 1
        # if counter == 10:
        #     break
    file.close()

    top_snp = {}
    top_snp['0.0-0.1'] = 0
    top_snp['0.1-0.2'] = 0
    top_snp['0.2-0.3'] = 0
    top_snp['0.3-0.4'] = 0
    top_snp['0.4-0.5'] = 0
    for snp in snp_dict:
        if float(snp_dict.get(snp)) - 0.5 <= 0.1:
            top_snp['0.0-0.1'] += 1
        if float(snp_dict.get(snp)) - 0.5 <= 0.2 and float(snp_dict.get(snp)) - 0.5 > 0.1:
            top_snp['0.1-0.2'] += 1
        if float(snp_dict.get(snp)) - 0.5 <= 0.3 and float(snp_dict.get(snp)) - 0.5 > 0.2:
            top_snp['0.2-0.3'] += 1
        if float(snp_dict.get(snp)) - 0.5 <= 0.4 and float(snp_dict.get(snp)) - 0.5 > 0.3:
            top_snp['0.3-0.4'] += 1
        if float(snp_dict.get(snp)) - 0.5 <= 0.5 and float(snp_dict.get(snp)) - 0.5 > 0.4:
            top_snp['0.4-0.5'] += 1
    labels = top_snp.keys()
    sizes = top_snp.values()
    explode = (0.05, 0.05, 0.05, 0.05, 0.05)
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
    plt.title('SNPs frequency difference from the 0.5\n' + os.path.basename(input_path)[11:-4])
    plt.show()
    print(top_snp)


def stacked_bar(input_path):
    snp_dicts = {}
    top_snps = {}
    for path in input_path:
        # range = os.path.basename(path)[11:-15]     # ALL_CPG
        range = os.path.basename(path)[11:-8]
        file = open(path)
        snp_dict = {}
        counter = 0
        for line in file:
            split_line = line.split(",")
            snp_dict[split_line[0]] = float(split_line[1].split("\n")[0])
            counter += 1
            # if counter == 10:
            #     break
        snp_dicts[range] = snp_dict
        file.close()

        top_snp = {}
        top_snp['40-50%'] = 0
        top_snp['30-40%'] = 0
        top_snp['20-30%'] = 0
        top_snp['10-20%'] = 0
        top_snp['0-10%'] = 0
        for snp in snp_dicts[range]:
            if float(snp_dicts[range].get(snp)) - 0.5 <= 0.1:
                top_snp['40-50%'] += 1
            if float(snp_dicts[range].get(snp)) - 0.5 <= 0.2 and float(snp_dicts[range].get(snp)) - 0.5 > 0.1:
                top_snp['30-40%'] += 1
            if float(snp_dicts[range].get(snp)) - 0.5 <= 0.3 and float(snp_dicts[range].get(snp)) - 0.5 > 0.2:
                top_snp['20-30%'] += 1
            if float(snp_dicts[range].get(snp)) - 0.5 <= 0.4 and float(snp_dicts[range].get(snp)) - 0.5 > 0.3:
                top_snp['10-20%'] += 1
            if float(snp_dicts[range].get(snp)) - 0.5 <= 0.5 and float(snp_dicts[range].get(snp)) - 0.5 > 0.4:
                top_snp['0-10%'] += 1
        # top_snp['Missed'] = 482421 - counter  # ALL_CPG
        top_snp['Missed'] = 523 - counter
        top_snps[range] = top_snp
    ind = np.arange(4)
    data = {}
    # 482421
    # ALL_CPG
    # data['40-50%'] = (top_snps['50']['40-50%'] / 482421 * 100.0, top_snps['100']['40-50%'] / 482421 * 100.0, top_snps['200']['40-50%'] / 482421 * 100.0, top_snps['500']['40-50%'] / 482421 * 100.0)
    # data['30-40%'] = (top_snps['50']['30-40%'] / 482421 * 100.0, top_snps['100']['30-40%'] / 482421 * 100.0, top_snps['200']['30-40%'] / 482421 * 100.0, top_snps['500']['30-40%'] / 482421 * 100.0)
    # data['20-30%'] = (top_snps['50']['20-30%'] / 482421 * 100.0, top_snps['100']['20-30%'] / 482421 * 100.0, top_snps['200']['20-30%'] / 482421 * 100.0, top_snps['500']['20-30%'] / 482421 * 100.0)
    # data['10-20%'] = (top_snps['50']['10-20%'] / 482421 * 100.0, top_snps['100']['10-20%'] / 482421 * 100.0, top_snps['200']['10-20%'] / 482421 * 100.0, top_snps['500']['10-20%'] / 482421 * 100.0)
    # data['0-10%'] = (top_snps['50']['0-10%'] / 482421 * 100.0, top_snps['100']['0-10%'] / 482421 * 100.0, top_snps['200']['0-10%'] / 482421 * 100.0, top_snps['500']['0-10%'] / 482421 * 100.0)
    # data['Missed'] = (top_snps['50']['Missed'] / 482421 * 100.0, top_snps['100']['Missed'] / 482421 * 100.0, top_snps['200']['Missed'] / 482421 * 100.0, top_snps['500']['Missed'] / 482421 * 100.0)

    data['40-50%'] = (top_snps['50']['40-50%'] / 523 * 100.0, top_snps['100']['40-50%'] / 523 * 100.0, top_snps['200']['40-50%'] / 523 * 100.0, top_snps['500']['40-50%'] / 523 * 100.0)
    data['30-40%'] = (top_snps['50']['30-40%'] / 523 * 100.0, top_snps['100']['30-40%'] / 523 * 100.0, top_snps['200']['30-40%'] / 523 * 100.0, top_snps['500']['30-40%'] / 523 * 100.0)
    data['20-30%'] = (top_snps['50']['20-30%'] / 523 * 100.0, top_snps['100']['20-30%'] / 523 * 100.0, top_snps['200']['20-30%'] / 523 * 100.0, top_snps['500']['20-30%'] / 523 * 100.0)
    data['10-20%'] = (top_snps['50']['10-20%'] / 523 * 100.0, top_snps['100']['10-20%'] / 523 * 100.0, top_snps['200']['10-20%'] / 523 * 100.0, top_snps['500']['10-20%'] / 523 * 100.0)
    data['0-10%'] = (top_snps['50']['0-10%'] / 523 * 100.0, top_snps['100']['0-10%'] / 523 * 100.0, top_snps['200']['0-10%'] / 523 * 100.0, top_snps['500']['0-10%'] / 523 * 100.0)
    data['Missed'] = (top_snps['50']['Missed'] / 523 * 100.0, top_snps['100']['Missed'] / 523 * 100.0, top_snps['200']['Missed'] / 523 * 100.0, top_snps['500']['Missed'] / 523 * 100.0)

    # colors = plt.cm.get_cmap('viridis', 5).colors
    colors = ['grey', 'grey','grey','grey']

    p1 = plt.bar(ind, data['40-50%'], 0.8, color=colors[0])
    p2 = plt.bar(ind, data['30-40%'], 0.8, bottom=np.array(data['40-50%']), color=colors[1])
    p3 = plt.bar(ind, data['20-30%'], 0.8, bottom=np.array(data['40-50%']) + np.array(data['30-40%']), color=colors[2])
    p4 = plt.bar(ind, data['10-20%'], 0.8, bottom=np.array(data['40-50%']) + np.array(data['30-40%']) + np.array(data['20-30%']), color=colors[3])
    p5 = plt.bar(ind, data['0-10%'], 0.8, bottom=np.array(data['40-50%']) + np.array(data['30-40%']) + np.array(data['20-30%']) + np.array(data['10-20%']), color=colors[4])
    # p6 = plt.bar(ind, data['Missed'], 0.35, bottom=np.array(data['40-50%']) + np.array(data['30-40%']) + np.array(data['20-30%']) + np.array(data['10-20%']))

    plt.grid(axis='y')

    #todo to understand, why values are not sumed
    # {'50': {'40-50%': 7, '30-40%': 9, '20-30%': 9, '10-20%': 15, '0-10%': 0, 'Missed': 483}, '100': {'40-50%': 12, '30-40%': 14, '20-30%': 18, '10-20%': 19, '0-10%': 0, 'Missed': 460},
# plt.title('SNPs frequency difference from the 0.5\n' + os.path.basename(input_path)[11:-4])
    plt.ylabel('%CpG')
    plt.xlabel('Range from CpG (bp)')
    plt.title('Probability to find different SNPs')
    plt.xticks(ind, top_snps.keys())
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), top_snps['50'].keys(), title='SNP minor allele frequency')
    plt.show()
    print(top_snps)

if __name__ == '__main__':
    input_path = sys.argv[1]
    create_diagram(input_path)
    # create_pie(input_path)
    # stacked_bar(sys.argv[1:])
    # /cs/cbio/ilia/project/SNP/SNPquality_50_0.5.csv /cs/cbio/ilia/project/SNP/SNPquality_100_0.5.csv /cs/cbio/ilia/project/SNP/SNPquality_200_0.5.csv /cs/cbio/ilia/project/SNP/SNPquality_500_0.5.csv
