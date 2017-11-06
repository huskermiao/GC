#!/usr/lib/python
#-*- utf-8 -*-

import vcf
from optparse import OptionParser

msg_usage = 'usage: %prog [-I] vcf [-O] map file name'
descr ='''transfer vcf file into map file used in GenotypeCorrector.'''
optparser = OptionParser(usage = msg_usage, description = descr)
optparser.add_option('-I', '--input', dest = 'vcffile',
                      help = 'input your vcf file')
optparser.add_option('-O', '--output', dest = 'mapfile',
                      help = 'output name of map file')
options, args = optparser.parse_args()


parse_gt={'0|0':'A','0|1':'X','1|0':'X','1|1':'B','.|.':'-'}

def fb2map(fbfile, mapfile):
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
        locus = '%s-%s'%(chrom, pos)
        gtline.append(locus)
        for j in i.samples:
            gt = j['GT']
            real_gt = parse_gt[gt]
            gtline.append(real_gt)
        line = '\t'.join(gtline)+'\n'
        mapfile.write(line)
    mapfile.close()
    vcffile.close()

if __name__ == "__main__":
    I = options.vcffile
    O = options.mapfile
    fb2map(I, O)
