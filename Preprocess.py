# -*- coding: UTF-8 -*-

"""
Preprocess raw genotype dataset before start correction
"""
import os.path as op
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from GC.apps.base import ActionDispatcher, OptionParser
from scipy.stats import chisquare as chi
from collections import defaultdict

def main():
    actions = (
        ('CStest', 'Remove SNPs with segregation distortions by testing ghenotype ratio of two homozygous using chisquare test.'),
        ('FilterMissing', 'Remove SNPs or/and samples with extremely high missing rate'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def getChunk(fn):
    df0_chr = defaultdict(int)
    chr_order = []
    with open(fn) as f:
        f.readline()
        for i in f:
            j = i.split()[0].split('-')[0]
            df0_chr[j] += 1
            if j in chr_order:
                pass
            else:
                chr_order.append(j)
    if len(chr_order) != len(set(chr_order)):
        sys.exit('Please check your marker name and sort them by chr name.')
    return chr_order, df0_chr

def genRatioPlot(s1, s2, chrom, befaft, a='A', b='B'):
    cond = (s1>=1) & (s2>=1)
    ratio = s1[cond]/s2[cond]
    ratio = ratio[ratio<20]
    ax = ratio.plot(kind='hist', bins=40, grid=True, alpha=0.75, edgecolor='k')
    ax.set_title('Ratio of %s:%s %s filteration in %s'%(a, b, befaft, chrom))
    ax.set_xlabel('Ratio')
    ax.set_ylabel('Counts')
    plt.savefig('Ration_%s%s_%s_%s.pdf'%(a, b, chrom, befaft))
    plt.clf()

def genMarkerMissPlot(series, chrom, befaft):
    ax = series.plot(kind='hist', bins=40, grid=True, alpha=0.75, edgecolor='k')
    ax.set_title('Missing rate of SNPs in %s %s filteration'%(chrom, befaft))
    ax.set_xlabel('Missing rate')
    ax.set_ylabel('Counts')
    plt.savefig('MissingRate_SNPs_%s_%s.pdf'%(chrom, befaft))
    plt.clf()

def genSampleMissPlot(series, befaft):
    ax = series.plot(kind='hist', bins=40, grid=True, alpha=0.75, edgecolor='k')
    ax.set_title('Missing rate distribution of samples %s filteration'%befaft)
    ax.set_xlabel('Missing rate')
    ax.set_ylabel('Counts')
    plt.savefig('MissingRate_Samples_%s.pdf'%(befaft))
    plt.clf()

def FMR(chunk_df, chrm, ho1, ho2, he, m, c):
    '''
    Filter Missing by Row
    '''
    if not chunk_df.index.is_unique:
        sys.exit('marker name is not unique in %s.'%chrom)
    df_tmp = chunk_df.replace([ho1, ho2, he, m], [0,1,2,3])
    sample_num = df_tmp.shape[1]
    missing_snp = df_tmp.apply(lambda x: (x==3).sum()/sample_num, axis=1)
    genMarkerMissPlot(missing_snp, chrm, 'before')
    good_snp_rate = missing_snp <= float(c)
    genMarkerMissPlot(missing_snp[good_snp_rate], chrm, 'after')
    good_snp = chunk_df.loc[good_snp_rate, :]
    return good_snp

def FMC(chunk_df, ho1, ho2, he, m):
    """
    Fliter Missing by Column
    """
    if not chunk_df.index.is_unique:
        sys.exit('marker name is not unique in %s.'%chrom)
    df_tmp = chunk_df.replace([ho1, ho2, he, m], [0,1,2,3])
    missing_sample = df_tmp.apply(lambda x: (x==3).sum(), axis=0)
    return missing_sample

def FilterMissing(args):
    """
    %prog FilterMissing genotype_matrix_input genotype_matrix_output 
    
      genotype_matrix_input: specify the input genotype matrix file. 
        Check the manual to maker sure the format is correct:
        (https://github.com/freemao/Genotype-corrector/wiki/Genotype-Corrector).

      genotype_matrix_output: specify the output file name after filteration.
    """
    p = OptionParser(FilterMissing.__doc__)
    p.set_common_opts()
    p.add_option('--mode', default='1', choices=('1', '2', '3'), 
        help = 'specify the filteration mode. 1: only filtre SNP markers. 2: only filter samples. 3: filter both markers and samples. mode 1 is recommend cause when you perform missing filteration on markers, the missing rate of samples will also be decreased.')
    p.add_option('--cutoff_snp', default='0.45',
        help = "specify the missing rate cutoff for SNPs. SNP's missing rate higher than this value will be removed")
    p.add_option('--cutoff_sample', default='0.6',
        help = "specify the missing rate cutoff for samples. Sample's missing rate higher than this value will be removed")
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    matrix_in, matrix_out, = args

    chr_order, df0_chr = getChunk(matrix_in)
    reader = pd.read_csv(matrix_in, delim_whitespace=True, index_col=0,  iterator=True)
 
    if opts.mode == '1':
        Good_SNPs = []
        total_flt_snp, chr_sample_missing = 0, []
        for chrom in chr_order:
            chunk = df0_chr[chrom]
            df_tmp1 = reader.get_chunk(chunk)
            marker_num, sample_num = df_tmp1.shape
            print('handling %s including %s markers and %s samples...'%(chrom, marker_num, sample_num))
            good_snp = FMR(df_tmp1, chrom, opts.homo1, opts.homo2, opts.hetero, opts.missing, opts.cutoff_snp)
            Good_SNPs.append(good_snp)
            missing_sample = FMC(good_snp, opts.homo1, opts.homo2, opts.hetero, opts.missing) 
            chr_sample_missing.append(missing_sample)
            total_flt_snp += good_snp.shape[0]
        sample_missing_num = pd.concat(chr_sample_missing, axis=1)
        sample_missing_rate = sample_missing_num.apply(lambda x: sum(x)/total_flt_snp, axis=1)
        genSampleMissPlot(sample_missing_rate, 'after')
        df1 = pd.concat(Good_SNPs)
        df1.to_csv(matrix_out, sep='\t', index=True)
        print('%s markers left after fitering.'%df1.shape[0])
    
    elif opts.mode == '2':
        total_snp, chr_sample_missing = 0, []
        for chrom in chr_order:
            chunk = df0_chr[chrom]
            total_snp += chunk
            df_tmp1 = reader.get_chunk(chunk)
            marker_num, sample_num = df_tmp1.shape
            print("calculating sample's missing rate in %s..."%(chrom))
            missing_sample = FMC(df_tmp1, opts.homo1, opts.homo2, opts.hetero, opts.missing)
            chr_sample_missing.append(missing_sample)
        sample_missing_num = pd.concat(chr_sample_missing, axis=1)
        sample_missing_rate = sample_missing_num.apply(lambda x: sum(x)/total_snp, axis=1)
        genSampleMissPlot(sample_missing_rate, 'before')
        cond_missing = sample_missing_rate < float(opts.cutoff_sample)
        good_samples = sample_missing_rate.index[cond_missing]
        sample_missing_Grate = sample_missing_rate[cond_missing]
        if len(sample_missing_Grate.index)!=0:
            genSampleMissPlot(sample_missing_rate[cond_missing], 'after')
            reader = pd.read_csv(matrix_in, delim_whitespace=True, index_col=0,  iterator=True)
            Good_Samples = []
            for chrom in chr_order:
                chunk = df0_chr[chrom]
                df_tmp = reader.get_chunk(chunk)[good_samples]
                print('outputing good samples in %s...'%(chrom))
                Good_Samples.append(df_tmp)
            df1 = pd.concat(Good_Samples)
            df1.to_csv(matrix_out, sep='\t', index=True) 
            print('Total %s samples were removed.'%(sample_num-len(good_samples)))
        else:
            print("no sample's missing rate less than the cutoff %s."%opts.cutoff_sample)

    else:
        total_snp, total_flt_snp, chr_sample_missing = 0, 0, []
        Good_SNPs = []
        for chrom in chr_order:
            chunk = df0_chr[chrom]
            total_snp += chunk
            df_tmp1 = reader.get_chunk(chunk)
            marker_num, sample_num = df_tmp1.shape
            print("calculating missing rate in %s..."%(chrom))
            good_snp = FMR(df_tmp1, chrom, opts.homo1, opts.homo2, opts.hetero, opts.missing, opts.cutoff_snp)
            Good_SNPs.append(good_snp)
            missing_sample = FMC(good_snp, opts.homo1, opts.homo2, opts.hetero, opts.missing) 
            chr_sample_missing.append(missing_sample)
            total_flt_snp += good_snp.shape[0]
        sample_missing_num = pd.concat(chr_sample_missing, axis=1)
        sample_missing_rate = sample_missing_num.apply(lambda x: sum(x)/total_flt_snp, axis=1)
        genSampleMissPlot(sample_missing_rate, 'before')
        cond_missing = sample_missing_rate < float(opts.cutoff_sample)
        good_samples = sample_missing_rate.index[cond_missing]
        if len(good_samples)!=0:
            genSampleMissPlot(sample_missing_rate[cond_missing], 'after')
            Good_Samples = []
            for ck, chrom in zip(Good_SNPs, chr_order):
                df_tmp = ck[good_samples]
                print('outputing good samples in %s...'%(chrom))
                Good_Samples.append(df_tmp)
            df1 = pd.concat(Good_Samples)
            df1.to_csv(matrix_out, sep='\t', index=True) 
            removed_snp = total_snp-total_flt_snp
            removed_sample = sample_num-len(good_samples)
            print('Total %s samples and %s markers were removed.'%(removed_sample, removed_snp))
        else:
            print("no sample's missing rate less than the cutoff %s."%opts.cutoff_sample)
            

def CStest(args):
    """
    %prog CStest genotype_matrix_input genotype_matrix_output

      genotype_matrix_input: specify the input genotype matrix file. 
        Check the manual to make sure the format is correct:
        (https://github.com/freemao/Genotype-corrector/wiki/Genotype-Corrector).

      genotype_matrix_output: specify the output file name after filteration.

    For F2 populations, test A:B against the expected ration of 1:1.
    For RIL populations, in addition to testing A:B against 1:1, SNPs with heterozygous proportion exceeding 50% of total genotype calls will also be removed.
    For BCFn population, test A:B against the expected ration of 3:1.
    """
    p = OptionParser(CStest.__doc__)
    p.set_common_opts()
    p.add_option('--population', default = 'RIL', choices=('F2', 'RIL', 'BCFn'),
        help = 'specify the population stype.')
    p.add_option('--degree', default = '3', choices=('1','2','3','4','5'),
        help = """set the chi square test cutoff strigency.
1: most strigent(will get rid of more markers). 5: least strigent(get ride of less markers).
You can choose the proper degree level by checking how many markers left under different degree levels.""")
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    matrix_in, matrix_out, = args

    cutoff = float('0.%s1'%('0'*int(opts.degree)))
    print('chi_square test pvalue cutoff: %s'%cutoff)

    ratio_dict = {'F2':(1,1), 'RIL':(1,1), 'BCFn':(3,1)}
    r1, r2 = ratio_dict[opts.population]
    
    chr_order, df0_chr = getChunk(matrix_in) 
    reader = pd.read_csv(matrix_in, delim_whitespace=True, index_col=0,  iterator=True)
    Filtered = []
    for chrom in chr_order:
        chunk = df0_chr[chrom]
        print('handling %s including %s markers...'%(chrom, chunk))
        df_tmp1 = reader.get_chunk(chunk)
        if not df_tmp1.index.is_unique:
            sys.exit('marker name is not unique in %s.'%chrom)
        df_tmp1_num = df_tmp1.replace([opts.homo1, opts.homo2, opts.hetero, opts.missing], [0,1,2,3])
        counts = df_tmp1_num.apply(lambda x: x.value_counts(), axis=1)
        obs_homo1, obs_homo2 = counts[0], counts[1]
        genRatioPlot(obs_homo1, obs_homo2, chrom, 'before', opts.homo1, opts.homo2)
        cond_min5 = (obs_homo1>5) & (obs_homo2>5)
        observe_homo1, observe_homo2 = obs_homo1[cond_min5], obs_homo2[cond_min5]
        observe_total = observe_homo1 + observe_homo2
        if opts.population == 'BCFn':
            expect_homo1, expect_homo2 = observe_total * 0.75, observe_total * 0.25
        elif opts.population == 'F2':
            expect_homo1, expect_homo2 = observe_total * 0.5, observe_total * 0.5
        elif opts.population == 'RIL':
            expect_homo1, expect_homo2 = observe_total * 0.5, observe_total * 0.5
        chi_value, p_value = chi([observe_homo1, observe_homo2], [expect_homo1, expect_homo2])
        df_tmp2 = pd.DataFrame(dict(zip(['ob_1', 'ob_2', 'chi', 'p'], [observe_homo1, observe_homo2, chi_value, p_value])))
        cond_pvalue = df_tmp2['p'] > cutoff
        genRatioPlot(df_tmp2['ob_1'][cond_pvalue], df_tmp2['ob_2'][cond_pvalue], chrom, 'after', opts.homo1, opts.homo2)
        cond_final = df_tmp2.index[cond_pvalue]
        df_tmp3 = df_tmp1.loc[cond_final, :]
        Filtered.append(df_tmp3)
    df1 = pd.concat(Filtered)
    df1.to_csv(matrix_out, sep='\t', index=True)
    print('%s markers left after fitering.'%df1.shape[0])


if __name__ == "__main__":
    main()
