#!/usr/lib/python
import random

def addmiss(mapfile, rate, outputfile):
    f0 = open(mapfile)
    firstline = f0.readline()
    ll = len(firstline.split()[1:])
    seq = ''
    locus_info = []
    for i in f0:
        k = i.split()
        seq += ''.join(k[1:])
        locus_info.append(k[0])
    seqls = list(seq)
    n = len(seq)
    print 'total %s elements'%n
    change_n = int(n*float(rate))
    print 'you will substitute %s'%change_n
    idx_ls = range(n)
    change_idx = random.sample(idx_ls, change_n)
    for j in change_idx:
        seqls[j] = '-'
    myseq = ''.join(seqls)
    myseq_ls = str2ls(myseq, ll)
    print 'myseq_ls: %s'%myseq_ls
    f1 = open(outputfile,'w')
    f1.write(firstline)
    for p1,p2 in zip(locus_info, myseq_ls):
        newline = p1 +'\t'+'\t'.join(list(p2))+'\n'
        f1.write(newline)
    f0.close()
    f1.close()

def str2ls(string, len_line):
    myls = []
    total_len = len(string)
    test = total_len%len_line
    if test == 0:
        lines = total_len/len_line
        print 'total len: %s\tlines:%s\tlength of \
line:%s'%(total_len,lines,len_line)
        start = 0
        end = len_line
        for i in range(lines):
            myseq = string[start: end]
            myls.append(myseq)
            start += len_line
            end += len_line
    else:
        print 'total length can not exact division the length of line...'
    return myls

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4:
        addmiss(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print '''usage: python add_miss.py map_file miss_rate output_file
This script only used for testing the performance of GC under different missing data. 
 
map_file, the simulated map file without any error genotypes and missing data
miss_rate, how many missing data you want to introduce to your simulated genotype data. 
output_file, the output map file. 

example:
python add_miss.py original.map 0.5 missing_0.5.map
'''











