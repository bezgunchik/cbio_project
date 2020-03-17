# import os
#
# with open("/cs/cbio/ilia/blocks19.tsv", mode="r") as blocks:
#     blocks.readline()
#     c = 0
#     for line in blocks:
#         c += 1
#         if (c < 500):
#             continue
#         line_lst = line.split()
#         cmd = os.popen("python3 /cs/cbio/ilia/project/vcf_reader/SNPsearcher.py 19 " + line[0] + " " + line[1] + " 1000").read()
#         print(cmd)
#
#
#     print(c)


#line = [0.046,0.97732574,0.019141665,0.02008774,0.010657556]
#line = [0.3315,0.061074354,0.22135834,0.25317737,0.46881828]
line = [0.9464,0.048028205,0.9740917,0.9856352,0.9827282]

result = []
for v in line:
    sigma = 0
    for other in line:
        sigma += abs(v - other)
    result.append(sigma/(len(line) - 1))
print(result)
print(max(result))