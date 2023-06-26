# -*- coding: utf-8 -*-
## follow the script Fst.py
## python2.7

import argparse
import pandas
import numpy as np
import sys
import socket
import os
import subprocess
import gzip
import time
import math
from operator import itemgetter
import itertools
from itertools import chain
from concurrent.futures import ProcessPoolExecutor

##################################################################################
#######################                 group              #######################
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

##################################################################################
#######################               sampling             #######################
##################################################################################
def split_window(chromID,start,end,windowsize,stepsize):
    overlapsize = windowsize - stepsize

    length = end - start + 1
    bin_num = max(int(math.ceil((length - overlapsize)*1.0 / stepsize)),1)
    ex_len = bin_num * stepsize + overlapsize
    ex_start = int(max(start - (ex_len-length)/2.0, 1.0))
    ex_end = int(end + (ex_len-length)/2.0)

    region = pandas.DataFrame(columns=['chr','start','end'])
    region['chr'] = [chromID] * bin_num
    region['start'] = [ex_start + num*stepsize for num in range(bin_num)]
    region['end'] = region['start'] + windowsize - 1

    return region

def group_marker(block, chromID):
    blockdf = block[block['chr'] == chromID].copy(); bimtmp = bim[bim['chr'] == chromID].copy()
    blockdf['markers'] = blockdf.apply(lambda x: list(bimtmp[(bimtmp['pos']>=x['start']) & (bimtmp['pos']<=x['end'])].index), axis=1)

    return blockdf

def sample_site(way,bstime,size,outputprefix):
    print('Do sampling...')
    chromlist = list(bim['chr'].unique())
    
    ## length of each chromosome
    chromparm = pandas.DataFrame(columns=['chr','start','end'])
    chromparm['chr'] = chromlist
    chromparm['start'] = chromparm['chr'].apply(lambda x: bim[bim['chr']==x]['pos'].min())
    chromparm['end'] = chromparm['chr'].apply(lambda x: bim[bim['chr']==x]['pos'].max())

    ## cut each chromosome into windows
    block = pandas.concat(list(chromparm.apply(lambda x: split_window(x['chr'],x['start'],x['end'],size,size),axis=1)),ignore_index=True)
    block['index'] = block['chr'].apply(str) +":"+ block['start'].apply(str) +":"+ block['end'].apply(str)

    ## markers of each block
    with ProcessPoolExecutor(max_workers = len(chromlist)) as p:
        block = pandas.concat(list(p.map(group_marker, [block]*len(chromlist), chromlist)),ignore_index=True)

    blockdic = dict(zip(list(block['index']), list(block['markers'])))
    blockvalue = blockdic.values()

    ## resample
    subprocess.call('mkdir -p '+outputprefix+'_resample',shell=True)
    
    ## cann't use multiprocess
    if way == 'block':
        for bs in range(1,bstime+1):
            resample = list(np.random.choice(np.array(blockdic.keys()),size=len(blockdic),replace=True))
            getmarker = list(itemgetter(*resample)(blockdic))
            pandas.DataFrame({'site':list(chain(*getmarker))}).to_csv(outputprefix+'_resample/resample.site.'+str(bs)+'.list.gz',header=None,index=None,compression='gzip')
    elif way == 'snp':
        for bs in range(1,bstime+1):
            getmarker = [np.random.choice(x,size=1) for x in blockvalue if len(x)>0]
            pandas.DataFrame({'site':list(chain(*getmarker))}).to_csv(outputprefix+'_resample/resample.site.'+str(bs)+'.list.gz',header=None,index=None,compression='gzip')
    else:
        pass

    print('Done.\n')

    return block.shape[0]

##################################################################################
#######################                 FST                #######################
##################################################################################
def fst_one_pair_one_resample(siteseed,fstdata,sitepath,outputprefix):
    if sitepath == 'NO':
        sites = pandas.read_csv(outputprefix+'_resample/resample.site.'+str(siteseed)+'.list.gz',header=None)
    else:
        sites = pandas.read_csv(sitepath+'/resample.site.'+str(siteseed)+'.list.gz',header=None)
    sites = sites[sites[0].isin(list(set(fstdata.index) & set(sites[0])))]
    sites = list(sites[0])

    sample_numerator = fstdata.loc[sites,'numerator'].values
    sample_denominator = fstdata.loc[sites,'denominator'].values

    sample_fst = max(sample_numerator.sum() / sample_denominator.sum(), 0.0)

    return siteseed, sample_fst

def fst_one_pair_all_resample(pair,bstime,nprocess,freqs,sitepath,outputprefix):
    pop1, pop2 = pair
    print("Calculating FST for {} vs. {} ...".format(pop1,pop2))

    data1 = pandas.read_csv(freqs+'/'+pop1+'.frq.gz',sep='\s+',usecols=['SNP','A1','MAF','NCHROBS'],index_col='SNP',dtype={'A1':'category'})
    data2 = pandas.read_csv(freqs+'/'+pop2+'.frq.gz',sep='\s+',usecols=['SNP','A1','MAF','NCHROBS'],index_col='SNP',dtype={'A1':'category'})
    ## remove sites with high missing rate
    site2drop = list(set(data1[data1['NCHROBS']<=1].index) | set(data2[data2['NCHROBS']<=1].index))
    if len(site2drop) >0:
        if len(site2drop) == data1.shape[0]:
            boot_res = pandas.DataFrame([(s,np.nan) for s in range(1,bstime+1)],columns=['siteID',pop1+':'+pop2])
            boot_res.set_index('siteID',inplace=True)
            return boot_res
        else:
            data1.drop(site2drop,inplace=True)
            data2.drop(site2drop,inplace=True)
    else:
        pass

    ## allele calibration
    data2['ref'] = data1['A1']
    data2['MAF'] = abs((data2['A1']!=data2['ref']).values - data2['MAF'].values)
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
    ## merge, then sampling
    if len(site2drop) >0:
        merge = globals()['bim'].drop(site2drop).copy()
        merge['numerator'] = numerator; merge['denominator'] = denominator
        merge.drop(['chr','pos'],axis=1,inplace=True)
    else:
        merge = pandas.DataFrame({'numerator':numerator, 'denominator':denominator}, index=list(globals()['bim'].index))

    with ProcessPoolExecutor(max_workers = nprocess) as pool2:
        boot_res = list(pool2.map(fst_one_pair_one_resample, range(1,bstime+1), [merge]*bstime, [sitepath]*bstime, [outputprefix]*bstime))
    
    boot_res = pandas.DataFrame(boot_res,columns=['siteID',pop1+':'+pop2])
    boot_res.set_index('siteID',inplace=True)
    boot_res[pop1+':'+pop2] = boot_res[pop1+':'+pop2].astype('float32')

    return boot_res

def fst_all_pair_all_resample(pairlist,bstime,nprocess,freqs,sitepath,outputprefix):
    pairnum = len(pairlist)
    nprocess1 = min(pairnum, nprocess/2+1)
    nprocess2 = max(nprocess/nprocess1+1, 2)

    print("Totally {} pairs, working with {}*{} processes...".format(str(pairnum), str(nprocess1), str(nprocess2)))
    print("It takes about {} run(s).".format(str(int(math.ceil(pairnum*1.0 / nprocess1) * math.ceil(bstime*1.0 / nprocess2)))))

    with ProcessPoolExecutor(max_workers = nprocess1) as pool1:
        pair_res = list(pool1.map(fst_one_pair_all_resample, pairlist, [bstime]*pairnum, [nprocess2]*pairnum, [freqs]*pairnum, [sitepath]*pairnum, [outputprefix]*pairnum))
    pair_res = pandas.concat(pair_res, axis=1)

    print('Output global FST after bootstrapping...')
    pair_res.to_csv(outputprefix+'.Bootstrap_FST.txt.gz',sep='\t',compression='gzip',na_rep='NA')

    quantile = pair_res.quantile([0.025,0.975]).T
    quantile.to_csv(outputprefix+'.Bootstrap_FST_95CI.txt',sep='\t',na_rep='NA')

    return 'Done.'

def main():
    ## input
    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", type=str, required=True, \
                        help="/path/to/bplink/file/prefix. No duplicated snpID in the bim file. Same FID and IID in the fam file.")
    parser.add_argument("--info", type=str, required=True, \
                        help="individual info file, samples to be used, <IID> <population ID>, 2 columns, no header, space or tab separated.")
    parser.add_argument("--pair", type=str, required=False, default='NO', \
                        help="population pairs, 2 columns, <pop1> <pop2>, no header, space or tab separated.")
    parser.add_argument("--freqs", type=str, required=True, \
                        help="directory for frq.gz files, generated while running Fst.py.")
    parser.add_argument("--sites", type=str, required=False, default='NO', \
                        help="directory for prepared #boot resample.site.#n.list.gz files, #n ranges from 1 to #boot, (generated by this script if not provided).")
    parser.add_argument("--thread", type=int, required=False, default=25, \
                        help="number of processes, default: 25 processes.")
    parser.add_argument("--size", type=int, required=False, default=100000, \
                        help="Block size for sites sampling if sites not provided, default: 100KB.")
    parser.add_argument("--sampling", type=str, required=False, default='block', choices=['block','snp','NO'],\
                        help="NO: with site files prepared; block: sample moving blocks, snp: random sample one SNP for each block.")
    parser.add_argument("--boot", type=int, required=False, default=100, \
                        help="numbers of bootstrapping, default: 100.")
    parser.add_argument("--out", type=str, required=False, default='./out', \
                        help="prefix for output files")
    args = parser.parse_args()

    ## logging
    with open(args.out+'.logfile','w') as log:
        log.write('python {}\n'.format(sys.argv[0]))
        log.write('{}--bfile    {}\n'.format(' '*8, args.bfile))
        log.write('{}--info     {}\n'.format(' '*8, args.info))
        log.write('{}--pair     {}\n'.format(' '*8, args.pair))
        log.write('{}--freqs    {}\n'.format(' '*8, args.freqs))
        log.write('{}--sites    {}\n'.format(' '*8, args.sites))
        log.write('{}--thread   {}\n'.format(' '*8, str(args.thread)))
        log.write('{}--size     {}\n'.format(' '*8, str(args.size)))
        log.write('{}--sampling {}\n'.format(' '*8, args.sampling))
        log.write('{}--boot     {}\n'.format(' '*8, str(args.boot)))
        log.write('{}--out      {}\n\n'.format(' '*8, args.out))
        
        log.write('Hostname: '+socket.gethostname()+'\n')
        log.write('Working directory: '+os.getcwd()+'\n')
        log.write('Start time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## groups, pairs, and regions
    sampleinfo, groups, pairs = getgroup(args.info, args.pair)

    with open(args.out+'.logfile','a') as log:
        log.write("Totally {} groups, {} pairs included.\n".format(str(len(groups)),str(len(pairs))))
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## position info, and sampling
    pos = list(pandas.read_csv(args.freqs+'/'+groups[0]+'.frq.gz',sep='\s+',usecols=['SNP'])['SNP'])
    globals()['bim'] = pandas.read_csv(args.bfile+'.bim',sep='\t',header=None,names=['chr','ID','gdis','pos','a0','a1'],usecols=['chr','ID','pos'],index_col='ID',dtype={'chr':'category','pos':'int32'})
    globals()['bim'] = globals()['bim'].loc[pos]

    # check duplicated SNPs in the bim file, not allowed while sampling SNPs
    if len(set(globals()['bim'].index)) != globals()['bim'].shape[0]:
        print('Duplicated SNP ID found. Plz check the input bim file.')
        exit()
    else:
        pass

    if args.sites == 'NO':
        blocknum = sample_site(args.sampling, args.boot, args.size, args.out)
        
        with open(args.out+'.logfile','a') as log:
            log.write("Totally {} variants included.\n".format(str(bim.shape[0])))
            log.write("Grouped all variants into {} blocks\n".format(str(blocknum)))
            log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')
    else:
        pass

    ## FST
    fstcalculation = fst_all_pair_all_resample(pairs, args.boot, args.thread, args.freqs, args.sites, args.out)

    with open(args.out+'.logfile','a') as log:
        log.write("Done FST calculation.\n")
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')
        log.write("Type ls {}* to get all the output.\n".format(args.out))

    print('Done.')
    print("Type ls {}* to get all the output.\n".format(args.out))
    print('Have a Nice Day!')

if __name__ == '__main__':
    main()
