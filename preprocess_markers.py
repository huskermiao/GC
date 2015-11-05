#!/usr/lib/python
from  subprocess import call
import random
'''
compress those makers that are very close into a single marker which is present!
'''
def bin_makers(mapfile, min_len, output_file, heteros_gt):
    print '1 cycle...'
    f1 = open(mapfile)
    all_lines = f1.readlines()
    f2 = open(output_file, 'w')
    fir_line = all_lines[0]
    f2.write(fir_line)
    min_interval = int(min_len)
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
        if fir_chr == sec_chr:
            if sec_pos-fir_pos<=min_interval and (first[0] not in binnedPos) :
                new_locus_name = fir.split()[0]+'-'+str(sec_pos)
                new_gt = combine_gt(fir_gt, sec_gt, heteros_gt)
                bin_line = new_locus_name + '\t' + new_gt + '\n'
                gtMatrix.append(bin_line)
                binnedPos.append(sec_chr+'-'+str(sec_pos)) # add marker bined to the pos set
            elif sec_pos-fir_pos>min_interval and first[0] not in binnedPos :
                gtMatrix.append(fir)
            else:
                pass
        else:
            gtMatrix.append(fir)
    last_line = second_flag[-1]
    if last_line.split()[0] not in binnedPos:
        gtMatrix.append(last_line)
    remainN = len(gtMatrix)
    print '%s markers binned!'%len(binnedPos)

    cycle_n = 2
    while True:
        print '%s cycles'%cycle_n
        N, gtMatrix = cycle_bin(gtMatrix,min_interval,heteros_gt)
        print '%s markers binned in this cycle.'%(N)
        if N == 0:
            break
        cycle_n += 1
    f2.writelines(gtMatrix)


def cycle_bin(gtMatrix,min_interval, heteros_gt):
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
        if fir_chr == sec_chr:
            if sec_pos-fir_pos<=min_interval and (first[0] not in binned_Pos):
                new_locus_name = \
first[0]+'-'+'-'.join(second[0].split('-')[1:]) #(chr1-100-150) + (300-350)
                new_gt = combine_gt(fir_gt, sec_gt, heteros_gt)
                bin_line = new_locus_name + '\t' + new_gt + '\n'
                Matrix_gt.append(bin_line)
                binned_Pos.append(second[0]) # add marker bined to the pos list
            elif sec_pos-fir_pos>min_interval and (first[0])not in binned_Pos:
                Matrix_gt.append(fir)
            else:
                pass
        else:
            Matrix_gt.append(fir)
    last_line = second_flag[-1]
    if last_line.split()[0] not in binned_Pos:
        Matrix_gt.append(last_line)
    N2 = len(Matrix_gt)

    return N1-N2, Matrix_gt




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

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 5:
        bin_makers(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print 'usage:\npython bin_markers.py map_file min_interval output_file hetegt_letter'
