#!/usr/lib/python

def calhrate(map_file, ho_a, ho_b, he_x):
    f1 = open(map_file)
    fir_line = f1.readline()
    sm_info = fir_line.split()[1:]
    sm_seq = map(list,zip(*(i.split()[1:] for i in f1)))
    total_x, total_a, total_b = 0, 0, 0
    for i,j in zip(sm_info,sm_seq):
        print 'for sample: %s'%i
        t_n = 0
        t_x_n = 0
        t_a_n = 0
        t_b_n = 0
        for m in j:
            t_n += 1
            if m == he_x:
                t_x_n += 1
            if m == ho_a:
                t_a_n += 1
            if m == ho_b:
                t_b_n += 1
        total_x += t_x_n
        total_a += t_a_n
        total_b += t_b_n
        print 'X rate: %s'%(t_x_n/float(t_n))
        print 'A rate: %s'%(t_a_n/float(t_n))
        print 'B rate: %s'%(t_b_n/float(t_n))
    T = float(total_x+total_a+total_b)
    print 'For whole population:\nA: %s\nB: %s\nX: %s'%\
(total_a/T, total_b/T,total_x/T)

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 5:
        calhrate(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print 'usage:\npython cal_hrate.py mapfile ho_a ho_b he_x'
