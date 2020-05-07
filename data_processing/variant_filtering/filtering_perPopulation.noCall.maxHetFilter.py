# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman
# These filters should be carried out on the per-population VCF and will check for:

1. sites where there are any no-call genotypes (removed); this stringency can be adjusted
2. sites where all called genotypes are heterozygotes (removed) (this is no longer done across all populations at once)

usage: use python 2.7
python filtering_bespokeFiltersAndChecks.py [infile full path] [outfile full path] [error file full path]
"""
import sys
import gzip
import argparse


parser = argparse.ArgumentParser(description='Count the number of 0/0 sites that have have called genotypes in at least [your projection value  / 2 ] or more individuals. Note that easy SFS projection values are in *haploids* Not Diploids')
parser.add_argument("--vcf",required=True,help="path to vcf file")
parser.add_argument("--outfile",required=True,help="path to output file")
parser.add_argument("--errorfile",required=True,help="path to error file")
parser.add_argument("--maxNoCallFrac",required=False,help="set a maximum no call fraction. 1.0 indicates no filter; 0.2 allows up to 20% missing data.",default="")
parser.add_argument("--maxHetFilter",required=True,help="desired maximum het filter; 1.0 means that any site that has 100% het calls is excluded; 0.75 means that any site with >75% het calls is excluded. Note that this is at the full sample level. You can do more fine-scale filtering per population in EasySFS.",default="")

args = parser.parse_args()
filepath=args.vcf #input file
outname= args.outfile #  output file
errorname= args.errorfile # error file
maxNoCallFrac= args.maxNoCallFrac
maxHetFilter = args.maxHetFilter


################ write header lines to new vcf file #############
# this then prints out the rest of the header lines
def main_vcf_check(inputvcfilename,outfilename,errorfilename,maxHetFilter,maxNoCallFrac):
    # open your files:
    inVCF = gzip.open(inputvcfilename, 'r')
    outVCF = open(outfilename, 'w')
    errorVCF = open(errorfilename, 'w')
    # set up a counter of no-pass (nop) or pass (p) sites
    counter_nop=0 # count no pass
    counter_p=0 # count pass
    # set up the header of the new vcf:
    for line0 in inVCF:
        if line0.startswith('#'):
            outVCF.write(line0)
            continue
        ### For all other non-header lines:
        line=line0.strip().split('\t') # this splits line by tabs
        #CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	[indivudals]
        mygenoinfo=line[9:]
        # get genotype info
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        myHomRef=allCalls.count("0/0") + allCalls.count("0|0") # number of hom ref gts
        myHet=allCalls.count("0/1") + allCalls.count("1/0") + allCalls.count("0|1")+ allCalls.count("1|0") # number of het gts
        myHomAlt=allCalls.count("1/1") + allCalls.count("1|1") # num of hom alt gts
        myCalled=myHomRef+myHet+myHomAlt # num of called genotypes
        myNoCalled=allCalls.count("./.") # num of no called genotypes
        # do a series of checks
        # check if there are more than X% no-called genotypes
        if float(myNoCalled) > round(float(len(allCalls))* float(maxNoCallFrac)):
            errorVCF.write('# More than ' + str(maxNoCallFrac) + ' fraction of genotypes are  no-call at this site\n')
            errorVCF.write(line0)
            counter_nop+=1
        # check if all calls are heterozygous at the population level
        # adding a max het filter
        elif myHet !=0 and myHet >= myCalled*float(maxHetFilter):
            errorVCF.write("# found a site with >="+str(float(maxHetFilter)*100)+"% heterozygous genotypes\n")
            errorVCF.write(line0)
            counter_nop+=1
        else:
            outVCF.write(line0)
            counter_p+=1
    # write out final information
    errorVCF.write('## This many lines failed a filter and werent printed' + '\t' + str(counter_nop) + '\n')
    errorVCF.write('## This many lines PASSED and were printed' + '\t' + str(counter_p) + '\n')
    errorVCF.close()
    outVCF.close()
    inVCF.close()

#### run the function ##########
main_vcf_check(filepath,outname,errorname,maxHetFilter,maxNoCallFrac)
# exit
sys.exit()
