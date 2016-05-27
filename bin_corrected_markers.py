#!/usr/lib/python
#-*- coding:utf-8 -*-

from  subprocess import call
import random
from optparse import OptionParser

msg_usage = '''python %prog [options]
example:
python bin_corrected_markers.py -i input_file -m 2 -o output_file -l - -c X'''
descr = '''DESCRIPTION: compress those very similar continuous makers into a single
representative marker.'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-i', '--genotype-matrix', dest = 'matrix_filename',
                     help = '''The genotype matrix file is tab delimited text
containing each sample's marker in each position. For more
detail about genotype matrix file please see our wiki page:
(https://github.com/freemao/Genotype-corrector/wiki/Genotype-Corrector).''')
optparser.add_option('-m', '--min_mismatch', dest = 'minimum_mismatch',
                     help = "The minimum mismatch between two continuous \n\
                     markers for bin.")
optparser.add_option('-o', '--output', dest = 'output_filename',
                     help = 'Write the binned results to this file.')
optparser.add_option('-l', '--miss_letter', dest = 'missing_data_letter',
                     help = "missing date character in your genotype \n\
                     matrix file. Usually is '-' or '.'.")
optparser.add_option('-c', '--hetero_letter', dest ='heterozygous_genotype_letter',
                     help = "Heterozygous genotype letter in your genotype \n\
                     matrix file. Usually is 'h' or 'X'.")
options, args = optparser.parse_args()
def bin_makers(mapfile, min_miss_n, output_file, gt_miss, heteros_gt):
    print '1 cycle...'
    f1 = open(mapfile)
    all_lines = f1.readlines()
    f2 = open(output_file, 'w')
    fir_line = all_lines[0]
    f2.write(fir_line)
    min_miss_num = int(min_miss_n)
    first_flag = all_lines[1:-1]
    second_flag = all_lines[2:]
    gtMatrix = [] #all lines except first line
    binnedPos = [] #chr-pos which has been bined in second flag
    for fir, sec in zip(first_flag, second_flag):
        first = fir.split()
        second = sec.split()
        fir_gt = first[1:]
        sec_gt = second[1:]
        fir_chr = first[0].split('-')[0]
        sec_chr = second[0].split('-')[0]
        fir_pos = int(first[0].split('-')[-1])
        sec_pos = int(second[0].split('-')[1])
        miss_num = cal_n(fir_gt, sec_gt, gt_miss)
        if fir_chr == sec_chr:
            if first[0] not in binnedPos:
                if miss_num <= min_miss_num:
                    new_locus_name = fir.split()[0]+'-'+str(sec_pos)
                    new_gt = combine_gt(fir_gt, sec_gt, heteros_gt)
                    bin_line = new_locus_name + '\t' + new_gt + '\n'
                    gtMatrix.append(bin_line)
                    binnedPos.append(second[0]) # add marker bined to the pos set
                else: gtMatrix.append(fir)
            else: pass
        else:
            if first[0] not in binnedPos: gtMatrix.append(fir)
            else: pass

    last_line = second_flag[-1]
    if last_line.split()[0] not in binnedPos:
        gtMatrix.append(last_line)
    remainN = len(gtMatrix)
    print '%s markers binned!'%len(binnedPos)

    cycle_n = 2
    while True:
        print '%s cycles'%cycle_n
        N, gtMatrix = cycle_bin(gtMatrix,min_miss_num,gt_miss, heteros_gt)
        print '%s markers binned in this cycle.'%(N)
        if N == 0:
            break
        cycle_n += 1
    new_Matrix = gen_binned_names(gtMatrix)
    f2.writelines(new_Matrix)


def cycle_bin(gtMatrix,min_miss_num, gt_miss, heteros_gt):
    N1 = len(gtMatrix)
    first_flag = gtMatrix[0:-1]
    second_flag = gtMatrix[1:]
    Matrix_gt = [] #all lines except first line
    binned_Pos = [] #chr-pos which has been bined in second flag
    for fir, sec in zip(first_flag, second_flag):
        first = fir.split()
        second = sec.split()
        fir_gt = first[1:]
        sec_gt = second[1:]
        fir_chr = first[0].split('-')[0]
        sec_chr = second[0].split('-')[0]
        fir_pos = int(first[0].split('-')[-1])
        sec_pos = int(second[0].split('-')[1])
        miss_num = cal_n(fir_gt, sec_gt, gt_miss)
        if fir_chr == sec_chr:
            if first[0] not in binned_Pos:
                if miss_num <= min_miss_num:
                    new_locus_name = first[0]+'-'+'-'.join(second[0].split('-')[1:])
                    new_gt = combine_gt(fir_gt, sec_gt, heteros_gt)
                    bin_line = new_locus_name + '\t' + new_gt + '\n'
                    Matrix_gt.append(bin_line)
                    binned_Pos.append(second[0]) # add marker bined to the pos list
                else: Matrix_gt.append(fir)
            else: pass
        else:
            if first[0] not in binned_Pos: Matrix_gt.append(fir)
            else: pass
    last_line = second_flag[-1]
    if last_line.split()[0] not in binned_Pos: Matrix_gt.append(last_line)
    N2 = len(Matrix_gt)
    return N1-N2, Matrix_gt

def cal_n(fir, sec, gt_miss):
    n = 0
    for i,j in zip(fir, sec):
        if i != gt_miss and j != gt_miss and i != j:
            n += 1
    return n


def combine_gt(gt_list1,gt_list2, heteros_gt):
    new_gt_ls = []
    for i,j in zip(gt_list1,gt_list2):
        if i == j:new_gt_ls.append(i)
        elif i == '-':new_gt_ls.append(j)
        elif j == '-':new_gt_ls.append(i)
        elif i == heteros_gt: new_gt_ls.append(i)
        elif j == heteros_gt: new_gt_ls.append(j)
        else:
            new_gt_ls.append(random.choice([i,j]))
    return '\t'.join(new_gt_ls)

def gen_binned_names(gtMatrix):
    f = open('bin.info', 'w')
    new_Matrix = []
    for i in gtMatrix:
        j = i.split()[0].split('-')
        if len(j) > 2:
            markers_num = len(j)-1
            GTs = '\t'.join(i.split()[1:])
            newline = '%s-%s(%s)\t%s\n'%(j[0],j[1],markers_num,GTs)
            new_Matrix.append(newline)
            f.write('%s-%s(%s): %s\n'%(j[0],j[1],markers_num,i.split()[0]))
        else:new_Matrix.append(i)
    return new_Matrix


if __name__ == "__main__":
    I = options.matrix_filename
    M = options.minimum_mismatch
    O = options.output_filename
    L = options.missing_data_letter
    C = options.heterozygous_genotype_letter
    if I and M and O and L and C:
        bin_makers(I,M,O,L,C)
    else:
        print 'Add -h to show help.'
