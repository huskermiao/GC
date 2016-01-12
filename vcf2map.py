#!/usr/lib/python
#-*- utf-8 -*-

import vcf
from optparse import OptionParser

msg_usage = 'usage: %prog [-M] max %missing data(uninclude) [-I] vcf \
[-O] map file name'
descr ='''transfer vcf file generated from freebayes into map file which used to construct
genetic map.'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-M', '--min_missing', dest = 'min_missing_data', \
default = 0.5,
                      help = 'Maximal missing data rate of each line in \
vcf.default is 0.5, if mr = 0.5, this line will be retained.')
optparser.add_option('-I', '--input', dest = 'vcffile',
                      help = 'Please input your vcf file from freebayes')
optparser.add_option('-O', '--output', dest = 'mapfile',
                      help = 'The minimal depth of each SNP sites')
options, args = optparser.parse_args()

chrs = []
def parse_locus(chr, pos):
    if chr not in chrs:
        chrs.append(chr)
        print chr
    return '%s-p%s'%(chr, pos)

def parse_gt(genotype, depth):
    if genotype == '0/0' and depth >= 3:
        return 'A'
    elif genotype == '0/1' and depth >= 5:
        return 'X'
    elif genotype == '1/1' and depth >= 3:
        return 'B'
    else:
        return '-'

def judgeline(gtline, num_samp, min_missing_rate):
        dash_rate = gtline.count('-')/float(num_samp)
        if dash_rate > min_missing_rate:
            return '-'
        if dash_rate <= min_missing_rate:
            return '\t'.join(gtline)+'\n'

def fb2map(fbfile, mapfile, min_missing_rate):
    '''min_dp should be int and min_missing_rate should be a float'''
    vcffile = open(fbfile)
    myvcf = vcf.Reader(vcffile)
    samples = myvcf.samples # samples list in this vcf file
    sam_num = len(samples)
    print '%s samples: %s'%(sam_num, samples)
    mapfile = open(mapfile, 'w')
    firstline = 'locus\t' + '\t'.join(samples)+'\n'
    mapfile.write(firstline)
    contents = [] #will be writen to map file
    for i in myvcf:
        gtline = []
        chrom = i.CHROM
        pos = str(i.POS)
        locus = parse_locus(chrom, pos)
        gtline.append(locus)
        for j in i.samples:
            gt = j['GT']
            try: dp = j['DP']
            except:
                print 'Sample %s. No DP value in %s %s, set as 0'%(j, chrom, pos)
                dp = 0
            real_gt = parse_gt(gt, dp)
            gtline.append(real_gt)
        lastline = judgeline(gtline, sam_num, min_missing_rate)
        if lastline == '-':
            pass
        else:
            mapfile.write(lastline)
    mapfile.close()
    vcffile.close()

if __name__ == "__main__":
    M = float(options.min_missing_data)
    I = options.vcffile
    O = options.mapfile
    fb2map(I, O, M)
