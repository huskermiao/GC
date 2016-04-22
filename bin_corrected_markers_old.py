#!/usr/lib/python

f2 = open('tmp.map', 'w')
def compress_markers(mapfile, diss_n, gt_miss, compressed_map, bin_info):
    f1 = open(mapfile)
    f3 = open(bin_info, 'w')
    f2.write(f1.readline())
    fir_marker = f1.readline()
    f2.write(fir_marker)
    p = fir_marker.split()
    fir_name = p[0]
    f3.write(fir_name)
    fir_chr = fir_name.split('-')[0]
    fir_pos = '-'.join(fir_name.split('-')[1:])
    fir_ls = p[1:]
    for i in f1:
        j = i.split()
        sec_name = j[0]
        sec_chr = sec_name.split('-')[0]
        sec_pos = '-'.join(sec_name.split('-')[1:])
        sec_ls = j[1:]
        if sec_chr == fir_chr:
            n = compare_markers(fir_ls, sec_ls, gt_miss)
            if n <= int(diss_n):
                f3.write('-'+sec_pos)
            if n > int(diss_n):
                f2.write(i)
                f3.write('\n'+sec_name)
                fir_marker = i
                p = fir_marker.split()
                fir_name = p[0]
                fir_chr = fir_name.split('-')[0]
                fir_pos = '-'.join(fir_name.split('-')[1:])
                fir_ls = p[1:]
        if sec_chr != fir_chr:
            f2.write(i)
            f3.write('\n'+sec_name)
            fir_marker = i
            p = fir_marker.split()
            fir_name = p[0]
            fir_chr = fir_name.split('-')[0]
            fir_pos = '-'.join(fir_name.split('-')[1:])
            fir_ls = p[1:]
    f1.close()
    f2.close()
    f3.close()
    conciseName(bin_info, compressed_map)

def compare_markers(fir, sec, gt_miss):
    n = 0
    for i,j in zip(fir, sec):
        if i != gt_miss and j != gt_miss and i != j:
            n += 1
    return n
def conciseName(bin_info, compressed_map):
    f1 = open('tmp.map')
    f2 = open(bin_info)
    f3 = open(compressed_map, 'w')
    f3.write(f1.readline())
    for i,j in zip(f1, f2):
        genos = '\t'.join(i.split()[1:])+'\n'
        s = j.rstrip().split('-')
        if len(s) > 2:
            name = s[0]+'-'+s[1]+'-'+s[-1]+'(%s)'%(len(s)-1)
        else:
            name = j.rstrip()
        newline = name + '\t' + genos
        f3.write(newline)
    f1.close()
    f2.close()
    f3.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 6:
        compress_markers(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print 'usage:\npython bin_corrected_markers.py mapfile diss_n miss_letter \
compressed_mapfile bin_info_file'



