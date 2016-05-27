#!/usr/lib/python
#-*- utf-8 -*-

from optparse import OptionParser

msg_usage = '''python %prog [options]
example:
python cleanup.py -i input_file -o output_file -a A -b B -x X -c -'''
descr = '''DESCRIPTION: Remove original genotypes with brackets.'''

optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-i', '--genotype-matrix', dest = 'matrix_filename',
                     help = '''corrected genotype matrix file containing\n\
                     origianl genotypes with brackets''')
optparser.add_option('-o', '--output', dest = 'output_filename',
                     help = 'output file name.')
optparser.add_option('-a', '--ref_gt', dest = 'ref_genotype',
                     help = "Homozygous 0/0 genotype letter in your genotype \n\
                     matrix file. Usually is 'a' or 'A'.")
optparser.add_option('-b', '--alt_gt', dest = 'alt_genotype',
                     help = "Homozygous 1/1 genotype letter in your genotype \n\
                     matrix file. Usually is 'b' or 'B'.")
optparser.add_option('-x', '--hetero_gt', dest = 'heterozygous_genotype',
                     help = "Heterozygous genotype letter in your genotype \n\
                     matrix file. Usually is 'h' or 'X'.")
optparser.add_option('-c', '--miss_gt', dest = 'missing_data_character',
                     help = "Missing data character in your genotype \n\
                     matrix file. Usually is '-' or '.'.")
options, args = optparser.parse_args()

def cleanup(correctedfile,output,A,B,X,C):
    f1 = open(correctedfile)
    allC = f1.read()
    all1 = allC.replace('(%s)'%A, '')
    all2 = all1.replace('(%s)'%B, '')
    all3 = all2.replace('(%s)'%X, '')
    all4 = all3.replace('(%s)'%C, '')
    f2 = open(output,'w')
    f2.write(all4)
    f1.close()
    f2.close()
if __name__ == "__main__":
    I = options.matrix_filename
    O = options.output_filename
    A = options.ref_genotype
    B = options.alt_genotype
    X = options.heterozygous_genotype
    C = options.missing_data_character
    if I and O and A and B and X and C:
        cleanup(I,O,A,B,X,C)
    else:
        print 'Add -h to show help.'


