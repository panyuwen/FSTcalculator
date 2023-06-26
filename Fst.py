# -*- coding: utf-8 -*-

## Note:
## python2.7
## plink1.9 is required, should be called by typing "plink1.9"
## bplink format with gender info included, required for the X chromosome

'''
Usage:
--bfile:  required, /bplink/file/prefix, 
                    No duplicated snpID in the bim file. 
                    Same FID and IID in the fam file, with/without gender info.
--info:   required, individual info file
                    only samples in the info file will be used, 
                    <IID> <population ID>, 2 columns, no header, space or tab separated.
--pair:   optional, population pairs, 
                    2 columns, <pop1> <pop2>, no header, space or tab separated. 
                    FST calculated only for population pairs in the pair file,
                    OR, pairwisely, according to the populations listed in the info file
--point:  optional, whether to output the FST value for each variant
                    choice: Y N
                    default: N, output only global FST values
--region: optional, calculate FST values for the given regions, rather than the whole data set
                    4 columns: <region ID> <chrom ID> <start physical pos> <end physical pos>
--thread: optional, number of processes,
                    default: 20
--out:    optional, prefix of output files,
                    default: ./out
'''


import argparse
import sys
import socket
import os
import time
import itertools
import numpy as np
import pandas
import subprocess
from concurrent.futures import ProcessPoolExecutor

##################################################################################
#######################            group & region          #######################
##################################################################################
def getgroup(infofile,pairfile):
    print('Processing sample infos...')    
    sampledf = pandas.read_csv(infofile,sep='\s+',header=None)
    sampledf[1] = sampledf[1].astype(str)

    if pairfile == 'NO':
        grouplist = list(sampledf[1].unique())
        pairlist = list(itertools.combinations(grouplist,2))
    else:
        pairdf = pandas.read_csv(pairfile,sep='\s+',header=None)
        pairdf[0] = pairdf[0].astype(str); pairdf[1] = pairdf[1].astype(str)
        grouplist = list(set(list(pairdf[0])+list(pairdf[1])))
        pairlist = [tuple(x) for x in list(pairdf.values)]

    print('Done.\n')
    return sampledf, grouplist, pairlist

def getregion(regionfile):
    print('Processing region infos...')
    if regionfile == 'one':
        regiondf = pandas.DataFrame()
    else:
        regiondf = pandas.read_csv(regionfile,sep='\s+',header=None)
    
    print('Done.\n')
    return regiondf

##################################################################################
#######################              frequency             #######################
##################################################################################
def freqcal(sampledf,grouplist,bfile,regiondf,outprefix):
    print('Calculating frequency(s)...')
    if os.path.isdir(outprefix+'_freqs'):
        pass
    else:
        os.mkdir(outprefix+'_freqs')

    count = []
    for group in grouplist:
        while ((not os.path.isfile(outprefix+'_freqs/'+group+'.info')) | (not os.path.isfile(outprefix+'_freqs/'+group+'.frq.gz'))):
            count.append(group)
            ## info file
            sampledf[sampledf[1]==group][[0,0,1]].to_csv(outprefix+'_freqs/'+group+'.info',sep=' ',header=None,index=None)
            ## script
            if regiondf.shape[0] == 0:
                subprocess.call('plink1.9 --bfile '+bfile+' --keep '+outprefix+'_freqs/'+group+'.info --freq gz --out '+outprefix+'_freqs/'+group,shell=True)
            else:
                subprocess.call('plink1.9 --bfile '+bfile+' --keep '+outprefix+'_freqs/'+group+'.info --extract range '+outprefix+'.region --freq gz --out '+outprefix+'_freqs/'+group,shell=True)
            time.sleep(1)
    
    print('Done.\n')
    return len(set(count))

##################################################################################
#######################      frequency differentiation     #######################
##################################################################################
def select_markers(fst_region_df, merge_df, chromID):
    fst_regions = fst_region_df.copy()
    fst_regions['tmp'] = fst_regions.apply(lambda x: merge_df[(merge_df['pos']>=x['start']) & (merge_df['pos']<=x['end'])], axis=1)

    return fst_regions

def fstcal_one_pair_all_region(pair,regiondf,point,outprefix,nprocess):
    pop1, pop2 = pair
    print("Calculating FST for {} vs. {} ...".format(pop1,pop2))
    
    data1 = pandas.read_csv(outprefix+'_freqs/'+pop1+'.frq.gz',sep='\s+',usecols=['SNP','A1','MAF','NCHROBS'],index_col='SNP',dtype={'A1':'category'})
    data2 = pandas.read_csv(outprefix+'_freqs/'+pop2+'.frq.gz',sep='\s+',usecols=['SNP','A1','MAF','NCHROBS'],index_col='SNP',dtype={'A1':'category'})
    ## remove sites with high missing rate
    site2drop = list(set(data1[data1['NCHROBS']<=1].index) | set(data2[data2['NCHROBS']<=1].index))
    if len(site2drop) >0:
        if len(site2drop) == data1.shape[0]:
            return pop1+":"+pop2, np.nan, np.nan, np.nan
        else:
            data1.drop(site2drop,inplace=True)
            data2.drop(site2drop,inplace=True)
    else:
        pass
    ## allele calibration
    data2['ref'] = data1['A1']
    data2['MAF'] = abs((data2['A1']!=data2['ref']).values - data2['MAF'].values)
    #data2['MAF'] = data2.apply(lambda x: x['MAF'] if x['A1']==x['ref'] else 1.0-x['MAF'],axis=1)
    ## frequency and sample size
    freq1 = data1['MAF'].values; num1 = data1['NCHROBS'].values
    freq2 = data2['MAF'].values; num2 = data2['NCHROBS'].values
    ## msg. msp, nc, fst
    aver = (freq1*num1 + freq2*num2) / (num1 + num2)
    msg = (num1*freq1*(1.0-freq1) + num2*freq2*(1.0-freq2)) * 1.0 / ((num1-1.0) + (num2-1.0))
    msp = (num1*(freq1-aver)**2 + num2*(freq2-aver)**2) / (2-1.0)
    nc = ((num1+num2) - (num1**2+num2**2) / (num1+num2)) / (2-1.0)

    numerator = msp - msg
    denominator = msp + (nc-1.0)*msg
    ## global Fst, for each region & also the genome-wide
    # genome-wide (for all included regions)
    fst_global = numerator.sum() / denominator.sum()
    fst_global = max(fst_global,0.0)

    ## regional and each SNP
    fst_region = pandas.DataFrame()
    fsts = pandas.Series()    
    
    if (regiondf.shape[0] == 0) & (point == 'N'):
        pass
    else:
        merge = globals()['bim'].copy() ## <chr> <pos>, index: snpID
        if len(site2drop) >0:
            merge.drop(site2drop,inplace=True)
        else:
            pass
        merge['numerator'] = numerator; merge['denominator'] = denominator
        
        ## regional
        if regiondf.shape[0] > 0:
            fst_region = regiondf.copy() ## <ID> <chr> <start> <end>
            chromlist = list(fst_region['chr'].unique())

            with ProcessPoolExecutor(max_workers = nprocess) as p:
                fst_region = pandas.concat(list(p.map(select_markers, [fst_region[fst_region['chr']==chromID] for chromID in chromlist], [merge[merge['chr']==str(chromID)] for chromID in chromlist], chromlist)))
            
            fst_region = fst_region['tmp'].apply(lambda x: max(x['numerator'].sum() / x['denominator'].sum(), 0.0))
        else:
            pass

        ## Fst for each SNP
        if point == 'Y':
            fsts = merge['numerator'] / merge['denominator']
            fsts[fsts<0] = 0.0
        else:
            pass

    return pop1+":"+pop2, fst_global, fst_region, fsts


def fstcal_all_pair_all_region(pairlist,regiondf,point,nprocess,outprefix):
    nprocess1 = min(nprocess, len(pairlist))
    nprocess2 = int(nprocess / nprocess1)

    print("Totally {} pairs, working with {} processes...".format(str(len(pairlist)),str(nprocess1)))
    print("It takes about {} run(s).".format(str(len(pairlist)/nprocess1)))

    with ProcessPoolExecutor(max_workers = nprocess1) as pool:
        res = list(pool.map(fstcal_one_pair_all_region, pairlist, [regiondf]*len(pairlist), [point]*len(pairlist), [outprefix]*len(pairlist), [nprocess2]*len(pairlist)))
    
    WorkingPairList= [unit[0] for unit in res]
    
    print('\nOutput pairwise global FST...')
    GlobalFST = pandas.DataFrame({'pair':WorkingPairList,'globalFST':[unit[1] for unit in res]})
    GlobalFST['pop1'] = GlobalFST['pair'].apply(lambda x: x.split(':')[0])
    GlobalFST['pop2'] = GlobalFST['pair'].apply(lambda x: x.split(':')[1])
    GlobalFST = GlobalFST[['pop1','pop2','globalFST']]
    GlobalFST.to_csv(outprefix+'.Global_pairwise_FST.txt',sep='\t',header=None,index=None,na_rep='NA')

    tmp = GlobalFST[['pop2','pop1','globalFST']].copy(); tmp.columns = ['pop1','pop2','globalFST']
    GlobalFST = pandas.concat([GlobalFST,tmp],ignore_index=True)
    GlobalFSTMat = GlobalFST.pivot(index='pop1',columns='pop2')['globalFST'].reset_index()
    GlobalFSTMat.columns.name = None
    GlobalFSTMat.index = GlobalFSTMat['pop1'].values; GlobalFSTMat.drop('pop1',axis=1,inplace=True)
    GlobalFSTMat.to_csv(outprefix+'.Global_matrix_FST.txt',sep='\t',na_rep='NA')

    if regiondf.shape[0] == 0:
        pass
    else:
        print('\nOutput regional pairwise global FST...')
        RegionFST = pandas.concat([regiondf] + [unit[2] for unit in res],axis=1)
        RegionFST.columns = ['ID','chr','start','end'] + WorkingPairList
        RegionFST.sort_values(by=['chr','start','end'],ascending=True,inplace=True)
        RegionFST.to_csv(outprefix+'.Region_Global_FST.txt.gz',sep='\t',index=None,compression='gzip',na_rep='NA')

    if point == 'Y':
        print('\nOuput pairwise FST values for all variants...\n')
        VariantFST = pandas.concat([globals()['bim']] + [unit[3] for unit in res],axis=1)
        VariantFST.columns = ['chr','pos'] + WorkingPairList
        VariantFST.sort_values(by=['chr','pos'],ascending=True,inplace=True)
        VariantFST.to_csv(outprefix+'.Variant_FST.txt.gz',sep='\t',compression='gzip',na_rep='NA')
    else:
        pass

    return 'Done'

##################################################################################
#######################                main                #######################
##################################################################################
def main():
    ## input
    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", type=str, required=True, \
                        help="/path/to/bplink/file/prefix, with gender info. No duplicated snpID in the bim file. Same FID and IID in the fam file.")
    parser.add_argument("--info", type=str, required=True, \
                        help="individual info file, samples to be used, <IID> <population ID>, 2 columns, no header, space or tab separated.")
    parser.add_argument("--pair", type=str, required=False, default='NO', \
                        help="population pairs, 2 columns, <pop1> <pop2>, no header, space or tab separated.")
    parser.add_argument("--point", type=str, required=False, default='N', choices=['Y','N'], \
                        help="whether to output the FST value for each variant, memory consuming if 'Y'(Yes).")
    parser.add_argument("--region", type=str, required=False, default='one', \
                        help="region file, 4 columns: <region ID> <chrom ID> <start physical pos> <end physical pos>, calculate FST for each region, default: regard the whole genome as one region.")
    parser.add_argument("--thread", type=int, required=False, default=20, \
                        help="number of processes, default: 20 processes.")
    parser.add_argument("--out", type=str, required=False, default='./out', \
                        help="prefix for output files")
    args = parser.parse_args()

    ## logging
    with open(args.out+'.logfile','w') as log:
        log.write('python {}\n'.format(sys.argv[0]))
        log.write('{}--bfile  {}\n'.format(' '*8, args.bfile))
        log.write('{}--info   {}\n'.format(' '*8, args.info))
        log.write('{}--pair   {}\n'.format(' '*8, args.pair))
        log.write('{}--point  {}\n'.format(' '*8, args.point))
        log.write('{}--region {}\n'.format(' '*8, args.region))
        log.write('{}--thread {}\n'.format(' '*8, str(args.thread)))
        log.write('{}--out    {}\n\n'.format(' '*8, args.out))
        
        log.write('Hostname: '+socket.gethostname()+'\n')
        log.write('Working directory: '+os.getcwd()+'\n')
        log.write('Start time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## groups and pairs, and regions
    sampleinfo, groups, pairs = getgroup(args.info,args.pair)
    pandas.DataFrame(pairs).to_csv(args.out+'_pairs.info',sep=' ',header=None,index=None)
    regions = getregion(args.region)
    if regions.shape[0] == 0:
        regioncount = 1
    else:
        regioncount = regions.shape[0]
        regions = regions[[0,1,2,3]]
        regions[[1,2,3,0]].to_csv(args.out+'.region',sep=' ',header=None,index=None)
        regions.columns = ['ID','chr','start','end']
        regions = regions.astype({'chr':'category'})
        regions.index = regions.apply(lambda x: "%s:%s:%s:%s" % (x['ID'],x['chr'],x['start'],x['end']),axis=1)

    with open(args.out+'.logfile','a') as log:
        log.write("Totally {} groups, {} pairs included.\n".format(str(len(groups)),str(len(pairs))))
        log.write("Analyzing {} region(s).\n".format(str(regioncount)))
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## frequency
    freqcount = freqcal(sampleinfo,groups,args.bfile,regions,args.out)
    
    with open(args.out+'.logfile','a') as log:
        log.write("Frequency of totally {} group(s) were calculated\n".format(str(freqcount)))
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    # position info
    if (regions.shape[0] == 0) & (args.point == 'N'):
        globals()['bim'] = pandas.DataFrame()
    else:
        print("Reading physical positions...")
        pos = list(pandas.read_csv(args.out+'_freqs/'+groups[0]+'.frq.gz',sep='\s+',usecols=['SNP'])['SNP'])
        globals()['bim'] = pandas.read_csv(args.bfile+'.bim',sep='\t',header=None,names=['chr','ID','gdis','pos','a0','a1'],usecols=['chr','ID','pos'],index_col='ID',dtype={'chr':'category','pos':'int32'})
        globals()['bim'] = globals()['bim'].loc[pos]

        with open(args.out+'.logfile','a') as log:
            log.write("Physical positions accessed.\n")
            log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')
        print("Done.\n")

    # FST
    fstcalculation = fstcal_all_pair_all_region(pairs,regions,args.point,args.thread,args.out)
    with open(args.out+'.logfile','a') as log:
        log.write("Done FST calculation.\n")
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')
        log.write("Type ls {}* to get all the outputs.\n".format(args.out))

    print('Done.')
    print("Type ls {}* to get all the output.\n".format(args.out))
    print('Have a Nice Day!')

if __name__ == '__main__':
    main()

