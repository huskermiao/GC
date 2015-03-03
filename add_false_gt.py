#!/usr/lib/python
import random

def addfaslegt(mapfile,a_rate,b_rate,h_rate, outputfile):
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
    idxls = range(n)
    print 'total %s elements'%n
    idxls_a = []
    idxls_b = []
    idxls_h = []
    for m,n in enumerate(seqls):
        if n == 'A': idxls_a.append(m)
        if n == 'B': idxls_b.append(m)
        if n == 'X': idxls_h.append(m)
    change_a = int(len(idxls_a)*float(a_rate))
    change_b = int(len(idxls_b)*float(b_rate))
    change_h = int(len(idxls_h)*float(h_rate))
    print 'you will substitute %s a'%change_a
    print 'you will substitute %s b'%change_b
    print 'you will substitute %s h'%change_h
    change_idx_a = random.sample(idxls_a, change_a)
    change_idx_b = random.sample(idxls_b, change_b)
    change_idx_h = random.sample(idxls_h, change_h)
    for j in change_idx_a:
        seqls[j] = 'B'
    for j in change_idx_b:
        seqls[j] = 'A'
    for j in change_idx_h:
        seqls[j] = 'AB'[random.randint(0,1)]
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
    if len(sys.argv) == 6:
        addfaslegt(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print 'usage: \npython add_false.py map_file a_rate b_rate h_rate output_file'
