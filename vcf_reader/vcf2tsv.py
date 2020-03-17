import sys
import os

if __name__ == '__main__':
    file_name = sys.argv[1]
    dir_name = os.path.dirname(file_name)
    # print(file_name[:-4])
    col_num = 6
    with open(file_name, mode='r') as vcf:
        for i in range (56):
            vcf.readline()
        with open(file_name[:-4] + '.tsv', 'w') as csvFile:
            # fields = vcf.readline().split()[:5]
            # writer = csv.DictWriter(csvFile, delimiter='\t', fieldnames=fields)
            # writer.writeheader()
            # for line in vcf:
            #     # writer.
            #     writer.write(line.split()[:5])
            # print("writing completed")

            # headers = vcf.readline()
            # print(vcf.readline().split(',')[:5])

            for line in vcf:
                splited = line.split()
                line_lst = splited[:5]
                freq_col = splited[-1].split(',')
                if (len(freq_col[-1]) == 1):
                    freq = freq_col[-2]
                else:
                    freq = freq_col[-1]
                line_lst.append(freq.split(';')[0])
                if len(line_lst[3]) == 1 and len(line_lst[4].split(',')[0]) == 1:
                    line_str = ""
                    for i in range (col_num - 1):
                        line_str += line_lst[i] + '\t'
                    line_str += line_lst[-1] + '\n'
                    csvFile.write(line_str)
            csvFile.close()
        vcf.close()
#TODO to understand a new file format
