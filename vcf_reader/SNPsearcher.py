import sys
import time



def main():
    """
    Searches SNP sites with in required region
    """
    t1 = time.time()
    chr_num = sys.argv[1]
    # start and stop not including bounds
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    range = int(sys.argv[4])
    delta_freq = float(sys.argv[5])
    searcher(chr_num, start, end, range, delta_freq)
    t3 = time.time()
    t3 = t3 - t1
    print("All time: ", t3)


def searcher(chr_num, start, end, range, delta_freq):
    """
    :param chr_num: Chromosome number
    :param start: Block start
    :param end: Block end
    :param range: Searching range from the bounds
    :param delta_freq: Frequency range from entropy maximum value
    :return:
    """
    with open("/cs/cbio/ilia/Documents/common_snp_all_20180418.tsv", mode='r') as snp_dict:
        line = snp_dict.readline().split('\t')
        while (line[0] != chr_num):
            line = snp_dict.readline().split('\t')
        while (int(line[1]) <= (start - range)):
            line = snp_dict.readline().split('\t')
        counter_pre = 0
        counter_post = 0
        while (int(line[1]) < start):
            if checkFreq(float(line[5]), delta_freq):
                print(line)
                counter_pre += 1
            line = snp_dict.readline().split('\t')
        while (int(line[1]) < end):
            line = snp_dict.readline().split('\t')
        while (int(line[1]) <= end + range):
            if checkFreq(float(line[5]), delta_freq):
                print(line)
                counter_post += 1
            line = snp_dict.readline().split('\t')
    snp_dict.close()
    print("Counter_pre:", counter_pre)
    print("Counter_post:", counter_post)


def checkFreq(freq, delta_freq):
    """
    Checks if given frequency is in a given range from entropy maximum
    :param freq: Frequency of given SNP
    :param delta_freq: Frequency range from entropy maximum value
    :return: True if given frequency is in a given range from entropy maximum;
              False otherwise
    """
    max_entropy = 0.5
    if abs(freq - max_entropy) <= delta_freq:
        return True
    return False


if __name__ == '__main__':
    main()