import sys
import os

#
# if __name__ == '__main__':
#     file_name = sys.argv[1]
#     dir_name = os.path.dirname(file_name)
#     # print(file_name[:-4])
#     col_num = 6
#     with open(file_name, mode='r') as vcf:
#         for i in range (56):
#             vcf.readline()
#         with open(file_name[:-4] + '.tsv', 'w') as csvFile:
#             # fields = vcf.readline().split()[:5]
#             # writer = csv.DictWriter(csvFile, delimiter='\t', fieldnames=fields)
#             # writer.writeheader()
#             # for line in vcf:
#             #     # writer.
#             #     writer.write(line.split()[:5])
#             # print("writing completed")
#
#             # headers = vcf.readline()
#             # print(vcf.readline().split(',')[:5])
#
#             for line in vcf:
#                 splited = line.split()
#                 line_lst = splited[:5]
#                 freq_col = splited[-1].split(',')
#                 if (len(freq_col[-1]) == 1):
#                     freq = freq_col[-2]
#                 else:
#                     freq = freq_col[-1]
#                 line_lst.append(freq.split(';')[0])
#                 if len(line_lst[3]) == 1 and len(line_lst[4].split(',')[0]) == 1:
#                     line_str = ""
#                     for i in range (col_num - 1):
#                         line_str += line_lst[i] + '\t'
#                     line_str += line_lst[-1] + '\n'
#                     csvFile.write(line_str)
#             csvFile.close()
#         vcf.close()
# #TODO to understand a new file format

PRIME_LOCATION = 5
OBSERVED_LOCATION = 6
with open("/cs/zbio/ilia/snp151_9_9_2019.ssv", mode='r') as vcf:
    with open("/cs/zbio/ilia/snp151_11_9_2019.tsv", mode='w') as tsvFile:
        for line in vcf:
            fields = line.split(" ")
            prime = fields[PRIME_LOCATION]
            observed = fields[OBSERVED_LOCATION].split("/")
            minor = "";
            if observed[0] == prime:
                minor = observed[1]
            else:
                minor = observed[0]

            primeFreq = "";
            if fields[7] > fields[8]:
                primeFreq = fields[7]
            else:
                primeFreq = fields[8]
            primeFreq = primeFreq.rstrip("\n")
            tsvFile.write(fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + fields[5] + "\t" + minor + "\t" + primeFreq + "\n")
            #print(fields[0] + " " + fields[1] + " " + fields[2] + " " + fields[3] + " " + fields[4] + " " + fields[5] + " " + minor + " " + primeFreq)
        tsvFile.close()
    vcf.close()

