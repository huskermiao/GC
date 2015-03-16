#!/usr/lib/python

def cal_correct_rate(ori_file, corrected_file):
    orif = open(ori_file)
    samples = orif.readline().split()[1:]
    tmp = samples.pop()
    print 'samples: %s'%samples
    o_gtseq_need_reverse = []
    total_markers = 0
    for i in orif:
        total_markers += 1
        j = i.rstrip().split()[1:]
        o_gtseq_need_reverse.append(j)
    o_gtseq = map(list, zip(*o_gtseq_need_reverse))
    original_seq = o_gtseq.pop()
    cf = open(corrected_file)
    cf.readline()
    c_gtseq_need_reverse = []
    for i in cf:
        j = i.rstrip().split()[1:]
        c_gtseq_need_reverse.append(j)
    c_gtseq = map(list, zip(*c_gtseq_need_reverse))
    tmp = c_gtseq.pop()
    for i, j, k in zip(samples, o_gtseq, c_gtseq):
        print 'tackling: %s'%i
        total_need_corrected = 0
        total_corrected = 0
        correct_r = 0
        correct_w = 0
        for m,n,o in zip(j,k,original_seq):
            if m != o:
                total_need_corrected += 1
            if m != n:
                total_corrected += 1
                if n[0] == o:
                    correct_r += 1
                if n[0] != o:
                    correct_w += 1
        print 'total_need_corrected: %s'%total_need_corrected
        print 'total_corrected: %s'%total_corrected
        print 'correct_r: %s'%correct_r
        print 'correct_w: %s'%correct_w
        error_rate = total_need_corrected/float(total_markers)
        correct_rate = total_corrected/float(total_need_corrected)
        correct_r_rate = correct_r/float(total_corrected)
        correct_w_rate = correct_w/float(total_corrected)
        print 'missing and error rate: %s'%error_rate
        print 'correct_rate: %s'%correct_rate
        print 'correct_right_rate: %s'%correct_r_rate
        print 'correct_false_rate: %s'%correct_w_rate
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        cal_correct_rate(sys.argv[1], sys.argv[2])
    else:
        print 'usage: cal_correct_rate.py original_map_file corrected_map_file'

