#!/usr/lib/python
import subprocess

def bin_makers(mapfile, min_len, output_file):
    child = subprocess.Popen('cp %s %s'%(mapfile,output_file),shell=True)
    child.wait()
    judge_over = [0]
    while True:
        lines = 0
        f1 = open(output_file)
        all_lines = f1.readlines()
        f2 = open(output_file, 'w')
        fir_line = all_lines[0]
        min_interval = int(min_len)
        f2.write(fir_line)
        lines += 1
        first_flag = all_lines[1:-1]
        second_flag = all_lines[2:]
        contents_ls = []
        pos_set = set()
        for fir, sec in zip(first_flag, second_flag):
            first = fir.split()
            second = sec.split()
            fir_gt = first[1:]
            sec_gt = second[1:]
            fir_chr = first[0].split('-')[0]
            sec_chr = second[0].split('-')[0]
            fir_pos = int(first[0].split('-')[-1])
            sec_pos = int(second[0].split('-')[-1])
            if fir_chr+'-'+str(fir_pos) not in pos_set:
                if (fir_chr == sec_chr and
                    sec_pos-fir_pos<min_interval and
                    judge_diff(fir_gt,sec_gt)):
                    new_locus_name = fir.split()[0]+'-'+str(sec_pos)
                    new_gt = combine_gt(fir_gt, sec_gt)
                    bin_line = new_locus_name + '\t' + new_gt + '\n'
                    f2.write(bin_line)
                    lines += 1
                    pos_set.add(sec_chr+'-'+str(sec_pos))
                else:
                    f2.write(fir)
                    lines += 1
            if fir_chr+'-'+str(fir_pos) in pos_set:
                continue
        last_line = second_flag[-1]
        if last_line.split()[0] not in pos_set:
            f2.write(last_line)
            lines += 1
        f1.close()
        f2.close()
        print lines
        if lines == judge_over[-1]:
            break
        if lines != judge_over[-1]:
            judge_over.append(lines)


def judge_diff(gt_list1, gt_list2):
    for i,j in zip(gt_list1, gt_list2):
        if i != '-' and j != '-' and i != j:
            return False
    return True

def combine_gt(gt_list1,gt_list2):
    new_gt_ls = []
    for i,j in zip(gt_list1,gt_list2):
        if i == j:new_gt_ls.append(i)
        elif i == '-' and j != '-':new_gt_ls.append(j)
        elif i != '-' and j == '-':new_gt_ls.append(i)
        else:print 'warning: different genotype !'
    return '\t'.join(new_gt_ls)

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4:
        bin_makers(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print 'usage:\npython bin_markers.py map_file min_interval output_file'
