#!/usr/lib/python

from scipy import stats
import random

def main(popu_n, chr_n, eachr_len, mark_dens, file_name, ho_a, ho_b, he_x):
    first_line = gen_1st_line(popu_n)
    locus_info = gen_markers_loci(chr_n, eachr_len, mark_dens)
    genotype_info = gen_F2_genotype(chr_n, eachr_len, mark_dens,\
popu_n,ho_a,ho_b,he_x)
    f0 = open(file_name, 'w')
    f0.write(first_line+'\n')
    need_trans_ls = []
    need_trans_ls.append(locus_info)
    for i in genotype_info:
        need_trans_ls.append(i)
    transed_ls = map(list, zip(*need_trans_ls))
    for l in transed_ls:
        new_l = '\t'.join(l)
        f0.write(new_l+'\n')
    f0.close()

def gen_1st_line(popu_n):
    SMs = '\t'.join([ 'SM-'+str(i+1) for i in range(popu_n)])
    firline = 'locus_name\t' + SMs
#    print 'firline:\n%s'%firline
    return firline

def gen_F2_genotype(chr_n, eachr_len, mark_dens, popu_n,ho_a,ho_b,he_x):
    F1_popu_1st_set_geno = []
    F1_popu_2nd_set_geno = []
    F2_popu = []
    for i in range(popu_n):
        first_set_geno, second_set_geno = \
gen_F1_genotype(chr_n,eachr_len,mark_dens,ho_a,ho_b)
        F1_popu_1st_set_geno.append(first_set_geno)
        F1_popu_2nd_set_geno.append(second_set_geno)
    index_sm = range(popu_n)
    for i in range(popu_n):
        pateranl = random.sample(F1_popu_1st_set_geno,1)
        materanl = random.sample(F1_popu_1st_set_geno,1)
#        print 'pateranl:\n%s'%pateranl
#        print 'materanl:\n%s'%materanl
        F2_sm_pa = []
        F2_sm_ma = []
        for pa,ma in zip(pateranl[0],materanl[0]):
            pa_num = stats.binom(1,0.5).rvs()
            ma_num = stats.binom(1,0.5).rvs()
            if pa_num == 1:
                pa_one = pa
            if pa_num == 0:
                pa_one = ''.join(change_ls(list(pa),ho_a, ho_b))
            if ma_num == 1:
                ma_one = ma
            if ma_num == 0:
                ma_one = ''.join(change_ls(list(pa),ho_a,ho_b))
            F2_sm_pa.append(pa_one)
            F2_sm_ma.append(ma_one)
        F2_sm_geno = []
        for i,j in zip(F2_sm_pa,F2_sm_ma):
            geno_seq = get_geno(i,j,ho_a,ho_b,he_x)
            F2_sm_geno.append(geno_seq)
        F2_popu.append(list(''.join(F2_sm_geno)))
#    print 'F2_popu:\n%s'%F2_popu
    return F2_popu

def get_geno(pa_one, ma_one, ho_a, ho_b, he_x):
    geno_ls = []
    for i,j in zip(pa_one, ma_one):
        if i == j == ho_a:geno_ls.append(ho_a)
        if i == j == ho_b:geno_ls.append(ho_b)
        if i != j:geno_ls.append(he_x)
    geno_seq = ''.join(geno_ls)
    return geno_seq

def gen_markers_loci(chr_n, eachr_len, mark_dens):
    '''eachr_len is a int list.
    [50, 45, 40, 56, 43, 35, 50, 45]
    '''
    chr_name = ['chr'+str(i+1) for i in range(chr_n)]
    each_chr_length = [int(i*1000000) for i in eachr_len]
    each_chr_marker_n = [int(i*mark_dens/1000) for i in each_chr_length]
    eachr_locus_ls = []
    for i,j,k in zip(each_chr_length, each_chr_marker_n, chr_name):
        pre_locus = random.sample(range(i),j)
        locus = [k+'-'+str(m) for m in sorted(pre_locus)]
        eachr_locus_ls.extend(locus)
#    print 'eachr_locus_ls:\n%s'%eachr_locus_ls
    return eachr_locus_ls

def gen_F1_genotype(chr_n, eachr_len, mark_dens, ho_a, ho_b):
    each_chr_length = [int(i*1000000) for i in eachr_len]
    each_chr_marker_n = [int(i*mark_dens/1000) for i in each_chr_length]
    first_set_geno = []
    second_set_geno = []
    for i in each_chr_marker_n:
        first_geno = ho_a*i
        second_geno = ho_b*i
        first_set_geno.append(first_geno)
        second_set_geno.append(second_geno)
#    print '1st_set_geno:\n%s'%first_set_geno
#    print '2nd_set_geno:\n%s'%second_set_geno
    exchange_n = stats.poisson.rvs(1,size=chr_n)
#    print 'exchange_n:\n%s'%exchange_n
    genoseq_idx = range(chr_n)
    for c1,c2,n,idx in zip(first_set_geno, second_set_geno, exchange_n, genoseq_idx):
        if n == 0:
#            print 'no exchange'
            continue
        else:
            c1_ls = list(c1)
            c2_ls = list(c2)
            geno_n = len(c1)
            middle_point = (geno_n-1)/2.0
            breakpoint_ls = random.sample(range(geno_n), n)
#            print 'breakpoint_ls:\n%s'%breakpoint_ls
            for i in breakpoint_ls:
                if i <= middle_point:
                    c1_ls[0:i+1] = change_ls(c1_ls,ho_a,ho_b)[0:i+1]
                    c2_ls[0:i+1] = change_ls(c2_ls,ho_a,ho_b)[0:i+1]
                else:
                    c1_ls[i:] = change_ls(c1_ls,ho_a,ho_b)[i:]
                    c2_ls[i:] = change_ls(c2_ls,ho_a,ho_b)[i:]
            new_c1 = ''.join(c1_ls)
            new_c2 = ''.join(c2_ls)
            first_set_geno[idx] = new_c1
            second_set_geno[idx] = new_c2
#    print '1st_set_geno:\n%s'%first_set_geno
#    print '2nd_set_geno:\n%s'%second_set_geno
    return first_set_geno, second_set_geno

def change_ls(the_ls, ho_a, ho_b):
    new_ls = []
    for i in the_ls:
        if i == ho_a:new_ls.append(ho_b)
        elif i == ho_b:new_ls.append(ho_a)
        else:
            print 'not a and also not b....'
    return new_ls

if __name__ == '__main__':
    import sys
    main(120, 6, [30,25,16,30,15,27,20,23], 0.04, 'F2_120i_6chr_0.03dens.map', 'A', 'B', 'X')
