#!/usr/lib/python

def stat(original,precorrect,corrected):
    f1 = open(original)
    f2 = open(precorrect)
    f3 = open(corrected)
    sm_line = f1.readline()
    f2.readline()
    f3.readline()
    samples = sm_line.split()[1:]
    print 'this population contain total % samples: \n%s'%(len(samples), samples)
    ori_ls = map(list, zip(*[i.split()[1:] for i in f1]))
    pre_ls = map(list, zip(*[i.split()[1:] for i in f2]))
    cor_ls = map(list, zip(*[i.split()[1:] for i in f3]))
    total_bases = len(samples)*len(ori_ls[0])
    print 'total bases: %s'%total_bases
    homo_base, hete_base, miss = 0, 0, 0
    for k in ori_ls:
        for l in k:
            if l == 'A' or l == 'B': homo_base += 1
            if l == 'X': hete_base += 1
            if l == '-': miss += 1
    print '\t%s (%s) homozygous genotypes'%(homo_base,homo_base/float(total_bases))
    print '\t%s (%s) heterozygous genotypes'%(hete_base,hete_base/float(total_bases))
    print '\t%s (%s) missing data'%(miss, miss/float(total_bases))
    homo_miss, homo_fals = 0, 0
    hete_miss, hete_fals = 0, 0
    for i, j in zip(ori_ls, pre_ls):
        for m,n in zip(i, j):
            if m != n:
                if (m == 'A' or m == 'B') and n == '-': homo_miss += 1
                if (m == 'A' or m == 'B') and n != '-': homo_fals += 1
                if m == 'X' and n == '-': hete_miss += 1
                if m == 'X' and n != '-': hete_fals += 1
    t_need_c = homo_miss+homo_fals+hete_miss+hete_fals
    print '%s genotypes needed to correct in this population'%(t_need_c)
    print '%s genotypes needed to correct in homozygous region'%(homo_miss+homo_fals)
    print '\t%s missing genotypes, %s wrong genotypes'%(homo_miss, homo_fals)
    print '%s genotypes needed to correct in heterozygous region'%(hete_miss+hete_fals)
    print '\t%s missing genotypes, %s wrong genotypes'%(hete_miss, hete_fals)

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
    t_c = cor_homo+cor_hete
    t_c_r = cor_r_homo+cor_r_hete
    t_c_w = cor_w_homo+cor_w_hete
    print '%s (correct rate: %s) genotypes were corrected, %s right (right rate:\
%s), %s wrong (wrong rate: %s)'\
%(t_c,t_c/float(t_need_c), t_c_r, t_c_r/float(t_c), t_c_w, t_c_w/float(t_c))
    print '\t%s genotypes were corrected in homo region, %s (%s)right, %s (%s)wrong'\
%(cor_homo,cor_r_homo,cor_r_homo/float(cor_homo),cor_w_homo,cor_w_homo/float(cor_homo))
    print '\t%s genotypes were corrected in hete region, %s (%s)right, %s (%s)wrong'\
%(cor_hete,cor_r_hete,cor_r_hete/float(cor_hete),cor_w_hete,cor_w_hete/float(cor_hete))

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4:
        stat(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print '''usage:\npython statistical_correct.py original.map precorrect.map corrected.map

The script only used for testing to compute the true positive rate, true negative rate, and accuracy. 

original.map: the original map file without any missing data and wrong genotypes
precorrect.map: the map file was artificially introduced the missing and wrong genotypes.
corrected.map: the corrected map file by GenotypeCorrector
'''



