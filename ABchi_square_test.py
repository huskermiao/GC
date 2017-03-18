from scipy.stats import chisquare
import numpy as np

def chitest(mapfile, gt_a, gt_b, p_value):
    p_value = float(p_value)
    expect_rate = np.array([0.5, 0.5])
    f = open(mapfile)
    firline = f.readline()
    print firline.rstrip()
    for i in f:
        a_count = i.count('\t'+gt_a)
        b_count = i.count('\t'+gt_b)
        observed = np.array([a_count, b_count])
        cvalue, pvalue = chisquare(observed, expect_rate*np.sum(observed))
        if pvalue > p_value:
            print i.rstrip()
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 5:
        chitest(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        print 'usage:\npython stats_check.py mapfile gt_a gt_b p_value'
