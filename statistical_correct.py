#!/usr/lib/python

def tongji(original,precorrect,corrected):
    f1 = open(original)
    f2 = open(precorrect)
    f3 = open(corrected)
    sm_line = f1.readline()
    f2.readline()
    f3.readline()
    samples = sm_line.split()[1:]
    print 'this population contain total % samples: %s'%(len(samples), samples)
    ori_ls = map(list, zip(*[i.split()[1:] for i in f1]))
    pre_ls = map(list, zip(*[i.split()[1:] for i in f2]))
    cor_ls = map(list, zip(*[i.split()[1:] for i in f3]))
    total_bases = len(samples)*len(ori_ls[0])
    print 'total %s bases'%total_bases
    homo_base, hete_base = 0, 0
    for k in ori_ls:
        for l in k:
            if l == 'A' or l == 'B': homo_base += 1
            if l == 'X': hete_base += 1
            if l == '-': print 'there are missing data in original file'
    print '\t%s homozygous genotype and %s heterozygous \
genotype'%(homo_base, hete_base)
    homo_miss, homo_fals = 0, 0
    hete_miss, hete_fals = 0, 0
    for i, j in zip(ori_ls, pre_ls):
        for m,n in zip(i, j):
            if m != n:
                if (m == 'A' or m == 'B') and n == '-': homo_miss += 1
                if (m == 'A' or m == 'B') and n != '-': homo_fals += 1
                if m == 'X' and n == '-': hete_miss += 1
                if m == 'X' and n != '-': hete_fals += 1
    print 'total %s genotypes needed to correct in this \
population'%(homo_miss+homo_fals+hete_miss+hete_fals)
    print '\tin homozygous region,there are %s genotypes needed to correct, \
contain %s missing genotypes, %s wrong genotypes'%(homo_miss+homo_fals, homo_miss, homo_fals)
    print '\tin heterozygous region,there are %s genotypes needed to correct, \
contain %s missing genotypes, %s wrong genotypes'%(hete_miss+hete_fals, hete_miss, hete_fals)

    cor_homo, cor_hete = 0, 0
    cor_w_homo, cor_w_hete = 0, 0
    cor_r_homo, cor_r_hete = 0, 0
    for i, j, k in zip(ori_ls, pre_ls, cor_ls):
        for m,n,o in zip(i, j, k):
            if n != o and m in 'AB' :
                cor_homo += 1
                if m == o[0]:cor_r_homo += 1
                else: cor_w_homo += 1
            if n != o and m == 'X' :
                cor_hete += 1
                if m == o[0]:cor_r_hete += 1
                else: cor_w_hete += 1
    print '%s genotypes were corrected, %s right, %s wrong'\
%(cor_homo+cor_hete,cor_r_homo+cor_r_hete,cor_w_homo+cor_w_hete)
    print 'in homo region, %s genotypes were corrected, %s right, %s wrong'\
%(cor_homo,cor_r_homo,cor_w_homo)
    print 'in hetero region, %s genotypes were corrected, %s right, %s wrong'\
%(cor_hete,cor_r_hete,cor_w_hete)







if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4:
        tongji(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print 'usage:\npython statistical_correct.py original precorrect \
corrected'
