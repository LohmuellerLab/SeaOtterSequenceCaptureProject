# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman

usage: python 2.7
python filtering_getNoCallPerInd.py [infile full path] [outfile full path]
"""
import sys
import gzip

filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file

################ this function will count the number of no-call genotypes per individual for a VCF ##############
def getNoCallPerInd(inputvcfilename,outfilename):
    # open your files:
    inVCF = gzip.open(inputvcfilename, 'r')
    outList = open(outfilename, 'w')
    # get sample names
    samples=[]
    for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break
    # reset vcf
    inVCF.seek(0)
    # set up an empty dictionary filled with zeroes to be a counter
    noCallDict=dict()
    for sample in samples:
        noCallDict[sample]=0
    # skip the header lines
    for line0 in inVCF:
        if line0.startswith('#'):
            continue
        ### For all other non-header lines:
        line=line0.strip().split('\t') # this splits line by tabs
        #CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	[indivudals]
        mygenoinfo=line[9:]
        # zip it together with samples
        # get genotype info
        ####### This is really useful: makes queriable dictionary
        # of sample names and their genotypes
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        callDict = dict(zip(samples,allCalls))
        # now want to detect any no calls and add to a counter for each individual
        # if a sample has a no call genotype, want to add 1 to its counter
        for sample in samples:
            if callDict[sample]=="./.":
                noCallDict[sample]+=1
    inVCF.close()
    outList.write('Sample\tNoCallCount\n')
    for sample, nocall in noCallDict.items():
        outList.write('{}\t{}\n'.format(sample, nocall))





#### run the function ##########
getNoCallPerInd(filepath,outname)

sys.exit()
