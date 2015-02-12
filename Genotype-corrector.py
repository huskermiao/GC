#!/usr/lib/python
#-*- coding:utf-8 -*-

from scipy import stats
from optparse import OptionParser
from skidmarks import wald_wolfowitz

msg_usage = 'usage: %prog [-I] input your map file [-O] output file name [-C] config file'
descr = '''Correct SNPs as genetic markers used for constructing genetic map...
Before correcting, remove those markers whose distance is less than the read
length, which ensure the random of h to a,b process.
Write a script to solve this problem.'''

optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-I', '--input', dest = 'inputfile',
                     help = 'Input your file need to tackle')
optparser.add_option('-O', '--output', dest = 'outputfile',
                     help = 'Your output file name')
optparser.add_option('-C', '--config', dest = 'configfile',
                     help = 'Your configuration file name')
options, args = optparser.parse_args()


def main(mapfile, outputfile, configfile):
    '''the mapfile can not contain any blank line
    win_size is int, >= will tackle, < not tackle'''

    po_type, gt_zeze, error_zeze, gt_zeon, error_zeon, gt_onon,\
error_onon, gt_miss, win_size = parseconfigfile(configfile)

    samples_ls, contigs_ls, indexes_ls, first_line = parse_mapfile_infos(mapfile)
    seqs_ls, loci_ls = parse_mapfile_seqs(mapfile, gt_zeze, gt_zeon, gt_onon)
    final_list_need_reverse = []
    for sam, seq in zip(samples_ls, seqs_ls):
        print '\nTackling sample: %s'%sam
#        print 'Its seq(including the genoytpe of all the contig):\n%s'%seq
        each_sm_seq_ls = []
        for ctg, idx in zip(contigs_ls, indexes_ls):
            print '\tTackling %s contig...'%ctg
            idx_st = int(idx.split('-')[0])
            idx_ed = int(idx.split('-')[1])
            orig_seq = seq[idx_st:idx_ed] #str
#            print '\tIts seq:\n%s'%orig_seq
#            h_seq, m_seq, t_seq = split_seq2hmt(orig_seq, win_size, gt_miss)
#            if len(m_seq) >= win_size:
            if len(orig_seq) >= win_size:
                orig_seq_no_h = h2a_b(orig_seq, gt_zeze, gt_zeon, gt_onon)
                windows_list = get_windows_list(orig_seq_no_h, win_size, gt_zeze,\
gt_zeon, gt_onon)
                scores_list = get_scores_list(windows_list, win_size, gt_miss,\
gt_zeze, gt_onon, gt_zeon)
                head_seq, main_seqlist, tail_seq = get_HTseq_Mseqlist(windows_list,\
scores_list, win_size, gt_zeze, gt_zeon, gt_onon, gt_miss, error_zeze,\
error_zeon, error_onon, po_type)
                main_seq = get_Mseq(main_seqlist, orig_seq_no_h, orig_seq, scores_list,gt_zeze, \
gt_zeon, gt_onon, gt_miss, win_size)
                final_seq = head_seq+main_seq+tail_seq

#            if len(m_seq) < win_size:
            if len(orig_seq) < win_size:
                final_seq = orig_seq
            final_seq_ls = compare_and_mark(orig_seq, final_seq)
            each_sm_seq_ls.extend(final_seq_ls)
        final_list_need_reverse.append(each_sm_seq_ls)
        write2file(final_list_need_reverse, outputfile, first_line, loci_ls)

def parse_mapfile_infos(mapfile):
    '''parse mapfile then return: samples list, contigs list, indexes list,
    first line'''
    f0 = open(mapfile)
    first_line = f0.readline()
    samples_list = first_line.split()[1:]
    contigs_list, indexes_list = [], []
    tmp_index = []
    for i, j in enumerate(f0):
        k = j.split()
#need change depend on the input format
        contigname = k[0].split('-')[0]
        if contigname not in contigs_list:
            tmp_index.append(i)
            contigs_list.append(contigname)
    tmp_index.append(i+1)
    f0.close()
    for i, j in zip(tmp_index[0:-1], tmp_index[1:]):
        indexes_list.append(str(i)+'-'+str(j))
    print 'There are total %s samples.'%len(samples_list)
    print 'samples_list:\n %s\n' %samples_list
    print 'There are total %s contigs.'%len(contigs_list)
    print 'contigs_list:\n %s\n' %contigs_list
    print 'Caution: later confirm the location of contig in the\
    name.Cp-AxS-S7-p160446 or s7-p0001 or ....'
#    print 'indexes_list:\n %s' %indexes_list
#    print 'first_line:\n %s' %first_line
    return samples_list, contigs_list, indexes_list, first_line

def parse_mapfile_seqs(mapfile, gt_zeze, gt_zeon, gt_onon):
    '''get each sample's sequence list and locus list which is first column of
    map file'''
    f0 = open(mapfile)
    firline = f0.readline()
    seqs_need_reverse = []
    loci_list = []
    for i in f0:
        seq = i.split()[1:]
        locus = i.split()[0]
        seqs_need_reverse.append(seq)
        loci_list.append(locus)
    f0.close()
    normal = map(list, zip(*seqs_need_reverse))
    seqs_list = []
    for i in normal:
        realseq = ''.join(i)
        seqs_list.append(realseq)
#    print 'seqs_list:\n %s'%seqs_list
#    print 'loci_list:\n %s'%loci_list
    return seqs_list, loci_list

'''
def split_seq2hmt(seq, win_size, gt_miss):
    head_base = seq[0]
    if head_base != gt_miss:
        h_seq = ''
        for i in seq:
            if i == head_base or i == gt_miss:
                h_seq += head_base
            else:
                break
    else:
        h_seq = ''

    tail_base = seq[-1]
    reverse_seq = seq[::-1]
    if tail_base != gt_miss:
        t_seq = ''
        for j in reverse_seq:
            if j == tail_base or j == gt_miss:
                t_seq += tail_base
            else:
                break
    else:
        t_seq = ''

    if len(h_seq) < win_size:
        pass
    if len(h_seq) >= win_size:
        h_seq = ''
    if len(t_seq) < win_size:
        pass
    if len(t_seq) >= win_size:
        t_seq = ''
    if len(t_seq) == 0:
        m_seq = seq[len(h_seq):]
    if len(t_seq) != 0:
        m_seq = seq[len(h_seq):][:-len(t_seq)]
#    print h_seq
#    print t_seq
#    print m_seq
    return h_seq, m_seq, t_seq
'''

def get_windows_list(m_seq, win_size, gt_zeze, gt_zeon, gt_onon):
    '''transfer the m_seq to windows list'''
    windows_list = []
    count_windows = len(m_seq)-win_size+1
    win_idx_st = 0
    win_idx_ed = win_size
    for i in range(count_windows):
        windows_list.append(m_seq[win_idx_st:win_idx_ed])
        win_idx_st += 1
        win_idx_ed += 1
#    print 'There are total %d windows in this contig orig_seq.'%count_windows
    print 'windows_list:\n %s'%windows_list
    return windows_list

def h2a_b(sequence, gt_zeze, gt_zeon, gt_onon):
    new_seq = ''
    for i in sequence:
        if i == gt_zeon:
            num = stats.binom(1, 0.5).rvs()
            if num == 1: new_base = gt_zeze
            if num == 0: new_base = gt_onon
        else:
            new_base = i
        new_seq += new_base
#    print 'old seq: %s'%sequence
#    print 'new seq: %s'%new_seq
    return new_seq

def get_scores_list(windows_list, win_size, gt_miss, gt_zeze, gt_onon, gt_zeon):
    '''calculate each window's score, return scores list
    if missing data ratio more than half of win_size, score is -.
    the score is represented by the ratio of a count and b count at each win'''
    expect_dash_count = (win_size+1)/2 #the bigest value of dash count
    a_count, b_count = 0, 0
    scores_list = []
    for win in windows_list:
#will change depend on map file format
        dash_count = win.count(gt_miss)
        if dash_count < expect_dash_count:
            for base in win:
#a, b, h, will also change depend on map file format
                if base == gt_zeze: a_count += 1
                if base == gt_onon: b_count += 1
            if b_count == 0:ratio = a_count/(b_count+0.1)
            if b_count != 0:ratio = a_count/float(b_count)
            a_count, b_count = 0, 0
        if dash_count >= expect_dash_count: ratio = '-'
        scores_list.append(ratio)
#    print 'scores_list:\n %s'%scores_list
    return scores_list

def get_HTseq_Mseqlist(windows_list, scores_list, win_size, gt_zeze,
gt_zeon,gt_onon,gt_miss,error_zeze, error_zeon, error_onon, po_type):
    '''H, T, M means head, tail, main'''
    middle_index = (win_size-1)/2
    miss_count_cutoff = (win_size+1)/2
    main_seq_ls = []
    for w, s in zip(windows_list, scores_list):
        b_count = w.count(gt_onon)
        ab_length = win_size - w.count(gt_miss)
#b_count, ab_length are used in binormal calculate, see the next function
        if s == '-': representive_base = w[middle_index]
        if s != '-': representive_base = get_expected_pro(b_count, ab_length,\
error_zeze, error_zeon, error_onon, po_type, gt_zeze, gt_zeon, gt_onon)
        main_seq_ls.append(representive_base) #list
    first_window = windows_list[0]
    last_window = windows_list[-1]
    if first_window.count(gt_miss) >= miss_count_cutoff:
        head_seq = first_window[0:middle_index]
    if last_window.count(gt_miss) >= miss_count_cutoff:
        tail_seq = last_window[middle_index+1:]
    if first_window.count(gt_miss) < miss_count_cutoff:
        head_seq = main_seq_ls[0]*middle_index
    if last_window.count(gt_miss) < miss_count_cutoff:
        tail_seq = main_seq_ls[-1]*middle_index
#    print 'main_seq:\n %s'%main_seq
#    print 'head_seq:\n %s'%head_seq
#    print 'tail_seq:\n %s'%tail_seq
    return head_seq, main_seq_ls, tail_seq

def get_expected_pro(b_count, ab_length, a_error, h_error, b_error, po_type,\
gt_zeze, gt_zeon, gt_onon):
    '''given the SNP error rates  of three genotypes, a_error, b_error, h_error,
    and the proportions of three genotypes in population.calculate the genotype
    which represent the sliding window'''
    if po_type == 'RIL':
#        print b_count,ab_length,a_error,h_error,b_error,po_type,gt_zeze,gt_zeon
        a_ex_prob = stats.binom.pmf(b_count, ab_length, a_error)*49.98
        h_ex_prob = stats.binom.pmf(b_count, ab_length, 0.5+h_error/2)*0.05
        b_ex_prob = stats.binom.pmf(b_count, ab_length, 1-b_error)*49.98
        if a_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_zeze
        if h_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_zeon
        if b_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_onon
#        print Genotype
        return Genotype
    elif po_type == 'F2':
        a_ex_prob = stats.binom.pmf(b_count, ab_length, a_error)*25
        h_ex_prob = stats.binom.pmf(b_count, ab_length, 0.5+h_error/2)*50
        b_ex_prob = stats.binom.pmf(b_count, ab_length, 1-b_error)*25
        if a_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_zeze
        if h_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_zeon
        if b_ex_prob == max(a_ex_prob, h_ex_prob, b_ex_prob):
            Genotype = gt_onon
#        print Genotype
        return Genotype
    else:
        print 'this tool just support RIL and F2 polulations now'

def get_Mseq(main_seqlist, orig_seq_no_h, orig_seq, scores_list, gt_zeze, gt_zeon, gt_onon,\
gt_miss, win_size):
    '''if win_size is 15, >=15 omit, <15 splited to a or b.
I suppose that the genotype of heterozygous is always right.
Maybe set this assumption as a argument?'''
    h_island_list, h_island_idx, h_island_score = get_h_islands_info(main_seqlist,\
    gt_zeon, gt_zeze, gt_onon, gt_miss,scores_list)
    middle_index = (win_size-1)/2
    for i in h_island_idx:
        j = [int(k) for k in i.split()]
        orig_start = j[0] + middle_index
        orig_end = j[-1] + middle_index + 1
        print 'orig_start: %s'%orig_start
        print 'orig_end: %s'%orig_end
        orig_seq_h = orig_seq[orig_start:orig_end]
        orig_seq_h2ab = orig_seq_no_h[orig_start:orig_end]
        print 'orig_seq_h2ab:%s'%orig_seq_h2ab
        print 'orig_seq_h:%s'%orig_seq_h
        print 'length of j: %s'%len(j)
        if len(j) <= middle_index:
            h_changed = []
            for idx in j:
                print 'j: %s'%j
                score = scores_list[idx]
                print 'score: %s'%score
                if score > 1:h_changed.append(gt_zeze)
                if score < 1:h_changed.append(gt_onon)
                print 'orig_seq_h: %s'%orig_seq_h
                if score == 1:h_changed.append(orig_seq[idx+middle_index])
                print 'h_changed: %s'%h_changed
            for p, q in zip(j ,h_changed):
                main_seqlist[p] = q
        if len(j) > middle_index:
            for_test_seq = ''.join(orig_seq_h2ab.split(gt_zeon)) #remove '-'
            ele_n = len(set(for_test_seq))
            if ele_n <= 1:
                h_changed = []
                for idx in j:
                    score = scores_list[idx]
                    if score > 1:h_changed.append(gt_zeze)
                    if score < 1:h_changed.append(gt_onon)
                    if score == 1:h_changed.append(orig_seq[idx+middle_index])
                for p, q in zip(j ,h_changed):
                    main_seqlist[p] = q
            if ele_n > 1:
                pro_p = wald_wolfowitz(for_test_seq)['p'] # > 0.05 or not
                if pro_p > 0.05: # this seq is random
                    print 'this seq is ramdom. h retained'
                if pro_p <= 0.05: # this seq is nonrandom
                    h_changed = []
                    for idx in j:
                        score = scores_list[idx]
                        if score > 1:h_changed.append(gt_zeze)
                        if score < 1:h_changed.append(gt_onon)
                        if score == 1:h_changed.append(orig_seq[idx+middle_index])
                    for p, q in zip(j ,h_changed):
                        main_seqlist[p] = q
    main_seq = ''.join(main_seqlist)
    print 'main_seq:\n%s'%main_seq
    return main_seq


def get_h_islands_info(main_seqlist,gt_zeon,gt_zeze,gt_onon,gt_miss,scores_list):
    tmp1_main_seqlist = main_seqlist[:]
#tmp1_main_seqlist = [h,h,a...b...b,h,h,h]
#    head_h_container = []
#head_h_container = [h,h]
#    for i in range(len(main_seqlist)):
#        if tmp1_main_seqlist[0] == gt_zeon:
#            head_h_container.append(tmp1_main_seqlist.pop(0))
#    tail_h_container = []
#tail_h_container = [h,h,h]
#    for i in range(len(main_seqlist)-len(head_h_container)):
#        if tmp1_main_seqlist[-1] == gt_zeon:
#            tail_h_container.append(tmp1_main_seqlist.pop())
#   for i in range(len(head_h_container)):
#        tmp1_main_seqlist.insert(0, gt_miss)
#    for i in range(len(tail_h_container)):
#        tmp1_main_seqlist.append(gt_miss)
#tmp1_main_seqlist = [-,-,a...b...b,-,-,-]
    tmp_seq1 = ''.join(tmp1_main_seqlist)
    print 'tmp_seq1: %s'%tmp_seq1
#tmp_seq1 = '--a...b...b---'
    h_indexlist = []
    for i, j in enumerate(tmp_seq1):
        if j == gt_zeon:
            h_indexlist.append(i)
    print 'h_indexlist: %s'%h_indexlist
    h_indexlist_next_use = h_indexlist[:]
#if tmp_seq1 = '--aaahhhhhbhhb--', h_indexlist = [5,6,7,8,9,11,12]
    h_island_list = tmp_seq1.replace(gt_zeze,' ').replace(gt_onon,' ')\
.replace(gt_miss,' ').split()
    print 'h_island_list: %s'%h_island_list
#h_island_list = [hhhhh, hh]'
    h_idx_island = []
    h_socre_island = []
    tmp_idx = ''
    tmp_score = ''
    for i in h_island_list:
        for j in i:
            tmp_idx += str(h_indexlist.pop(0))+'\t'
            tmp_score += str(scores_list[h_indexlist_next_use.pop(0)])+'\t'
        h_idx_island.append(tmp_idx)
        h_socre_island.append(tmp_score)
        tmp_idx = ''
        tmp_score = ''
    print 'h_idx_island: %s'%h_idx_island
#h_idx_island = ['0\t1\t2\t3\t4\t5\t', '10\t11\t']
    print 'h_socre_island: %s'%h_socre_island
#h_socre_island = [each h's score joined by \t]
    return h_island_list, h_idx_island, h_socre_island

def compare_and_mark(orig_seq, prefinal_seq):
    '''the genotype corrected wii add a star'''
    final_seq_ls = []
    for i, j in zip(orig_seq, prefinal_seq):
        if i == j:final_seq_ls.append(j)
        if i != j:final_seq_ls.append(j+'*')
#    print 'final_seq_ls:\n%s'%final_seq_ls
    return final_seq_ls

def write2file(need_reverse_ls, outputfile, first_line, loci_ls):
    reversed_ls = map(list, zip(*need_reverse_ls))
    f0 = open(outputfile, 'w')
    f0.write(first_line)
    for loc, gp in zip(loci_ls, reversed_ls):
        gpline = '\t'.join(gp)+'\n'
        f0.write(loc+'\t'+gpline)
    f0.close()

def parseconfigfile(configfile):
    f = open(configfile)
    namedict = {}
    for i in f:
        if i.startswith('#'):pass
        else:
            if i.split(): #judge is blank line or not
                j = i.strip().split(':')
                namedict[j[0].strip()]=j[1].strip()
    print 'Your parameters in configuration file:\n'
    print '\tPopulation type: %s'%namedict['Population_type']
    print '\tLetter for 0/0: %s'%namedict['Letter_for_0/0']
    print '\tLetter for 0/1: %s'%namedict['Letter_for_0/1']
    print '\tLetter for 1/1: %s'%namedict['Letter_for_1/1']
    print '\tCharacter for missing data: \
%s'%namedict['Character_for_missing_data']
    print '\tSliding window size: %s'%namedict['Sliding_window_size']
    po_type = namedict['Population_type']
    gt_zeze = namedict['Letter_for_0/0']
    error_zeze = namedict['error_rate_for_0/0']
    print '\tSNP error rate for 0/0: %s'%error_zeze
    gt_zeon = namedict['Letter_for_0/1']
    error_zeon = namedict['error_rate_for_0/1']
    print '\tSNP error rate for 0/1: %s'%error_zeon
    gt_onon = namedict['Letter_for_1/1']
    error_onon = namedict['error_rate_for_1/1']
    print '\tSNP error rate for 1/1: %s'%error_onon
    gt_miss = namedict['Character_for_missing_data']
    win_size = namedict['Sliding_window_size']
    return po_type, gt_zeze, float(error_zeze), gt_zeon, float(error_zeon), gt_onon,\
float(error_onon), gt_miss, int(win_size)

if __name__ == "__main__":
    I = options.inputfile
    O = options.outputfile
    C = options.configfile
    main(I, O, C)
