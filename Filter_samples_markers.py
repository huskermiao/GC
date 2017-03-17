#!/usr/lib/python
#-*- coding:utf-8 -*-

from optparse import OptionParser

msg_usage = '''python %prog [options]
example:
python Filter_samples_markers.py -i input -m 0.4 -s 0.5 -o output_file'''
descr = '''DESCRIPTION: Remove too bad samples(missing rate > 50%) and markers
(missing rate > 40%). The characters for 0/0,1/1,0/1, and missindg data are A,
B. X, and - respectly.
'''

optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-i', '--genotype-matrix', dest = 'matrix_filename',
                     help = '''The genotype matrix file is tab delimited text
containing each sample's marker in each position. For more detail about genotype
matrix file please see our wiki page:
(https://github.com/freemao/Genotype-corrector/wiki/Genotype-Corrector).''')
optparser.add_option('-m', '--marker_miss', dest = 'Marker_missing_rate',
                     help = "The cutoff of markers' missing rate.")
optparser.add_option('-s', '--sample_miss', dest = 'Sample_missing_rate',
                     help = "The cutoff of samples' missing rate.")
optparser.add_option('-o', '--output', dest = 'output_filename',
                     help = 'Write the filtered results to this file.')
options, args = optparser.parse_args()

def filterMandS(mapfile, marker_missrate, sample_missrate, output):
    f1 = open(mapfile)
    firstColumn = []
    needConverted = []
    for i in f1:
        j = i.split()
        firstColumn.append(j[0])
        needConverted.append(j[1:])
    converted = map(list, zip(*needConverted))

    filteredList = []
    filteredList.append(firstColumn)
    for i in converted:
        N = len(i)-1
        miss_n = i.count('-')
        sample_miss_rate = miss_n/float(N)
        if sample_miss_rate < float(sample_missrate):
            filteredList.append(i)
    converted = map(list, zip(*filteredList))

    f2 = open(output, 'w')
    f2.write('\t'.join(converted[0])+'\n')
    for i in converted[1:]:
        newline = '\t'.join(i)+'\n'
        miss_n = i.count('-')
        marker_miss_rate = miss_n/float(len(i)-1)
        if marker_miss_rate < float(marker_missrate):
            f2.write(newline)
    f1.close()
    f2.close()

if __name__ == "__main__":
    I = options.matrix_filename
    M = options.Marker_missing_rate
    S = options.Sample_missing_rate
    O = options.output_filename
    if I and M and S and O:
        filterMandS(I,M,S,O)
    else:
        print 'Add -h to show help.'
