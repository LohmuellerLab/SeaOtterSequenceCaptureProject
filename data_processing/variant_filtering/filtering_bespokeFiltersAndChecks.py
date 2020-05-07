# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman
# These filters will check for:
# 1. ref or alt alleles that aren't a single letter (AGCT) or . (alt)
# 2. genotypes that aren't in 0/0, 0/1, 1/1 or ./. (maybe it's phased, etc)
# 3. must have qual score
# 4. must be PASS for site
# 5. make sure not missing DP, AN, GT, AD, DP or GQ/RGQ
# 6. make sure no called genotype is missing any info from the genotype info field
# 7. get rid of sites where all calls are 0/1 (all hets) [20180901: note you must do this again per-population as well]
# X. *doesn't* get rid of sites with >20% called sites [20180901: because I moved this to a per-population script]
# 8. instead, gets rid of sites where >90% are uncalled (should have been done already in step 1a, but if it wasn't, do it here.)
# this is very lenient ON PURPOSE because don't want to be stringent across populations. Later on can be stringent.
# 9. updates AN and AC based on final sets of calls (these aren't updated when GATK does genotype filtering)

# this script does NOT: change any genotypes; do any genotype filtering; change any FT fields for genotypes (./. gts will still be PASS if they started as ./. -- bit of GATK weirdness that isn't fatal)

usage: use python 2.7
python filtering_bespokeFiltersAndChecks.py [infile full path] [outfile full path] [error file full path]
"""
import sys
import gzip

filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file
errorname= sys.argv[3] # error file
maxNoCallFrac=sys.argv[4] # this is fraction of no-call genotypes you'll permit per site (anything greater than that will be excluded)

################ write header lines to new vcf file #############
# this then prints out the rest of the header lines
def main_vcf_check(inputvcfilename,outfilename,errorfilename):
    # open your files:
    inVCF = gzip.open(inputvcfilename, 'r')
    outVCF = open(outfilename, 'w')
    errorVCF = open(errorfilename, 'w')
    # set up your possible ref and alt bases
    refbases=set(['A', 'C', 'G', 'T'])
    altbases=set(['A', 'C', 'G', 'T', '.'])
    # set up your expected genotypes (change if phased!!!!)
    expected_genos=set(['0/0', '0/1', '1/1', './.']) # this doesn't expect phase
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
        myref=line[3]
        myalt=line[4]
        myqual=line[5]
        myfilter=line[6]
        myinfo=line[7]
        infoFields=dict(s.split('=') for s in myinfo.split(";")) # this makes my info into a dictionary that you can query; e.g. infoFields['AN'] will give you the AN value
        myformat=line[8]
        mygenoinfo=line[9:]
        # get genotype info
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        myHomRef=allCalls.count("0/0") + allCalls.count("0|0") # number of hom ref gts
        myHet=allCalls.count("0/1") + allCalls.count("1/0") + allCalls.count("0|1")+ allCalls.count("1|0") # number of het gts
        myHomAlt=allCalls.count("1/1") + allCalls.count("1|1") # num of hom alt gts
        myCalled=myHomRef+myHet+myHomAlt # num of called genotypes
        myNoCalled=allCalls.count("./.") # num of no called genotypes
        # get AC and AN values:
        myAC=myHet+(2*myHomAlt) # alternate alleles
        myAN=myCalled *2 # all called alleles
        # do a series of checks
        # check if ref base is a single letter A, G C or T (gets rid of weird indel refs)
        if myref not in refbases:
        		errorVCF.write('# ref not AGCT\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        # check if alternative allele is a proper single letter or .
        elif myalt not in altbases:
        		errorVCF.write('# alt not . or AGCT\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        # make sure that all genotypes are either 0/0, 0/1, 1/1 or ./.
        elif True in [x not in expected_genos for x in allCalls]:
                errorVCF.write('# Genotypes are not in an expected format -- are they phased?\n')
                errorVCF.write(line0)
                counter_nop+=1
        # make sure there's a qual score
        elif myqual == '.':
        		errorVCF.write('# QUAL = .\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        # make sure there are only PASS sites
        elif myfilter != 'PASS':
        		errorVCF.write('# FILTER is not PASS\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        # Make sure not missing DP, AN, GT, AD, DP or GQ/RGQ
        elif 'DP' not in myinfo or 'AN' not in myinfo:
        		errorVCF.write('# site info has no DP or AN\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        elif 'GT' not in myformat or 'DP' not in myformat or 'GQ' not in myformat:
        		errorVCF.write('# GT, DP or GQ (which also does RGQ) not in FORMAT field\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        # changed AD filter; only care if a SNP is missing AD, not if non var are missing it
        elif 'AD' not in myformat and myAC!=0:
        		errorVCF.write('# AD Missing from FORMAT field\n')
        		errorVCF.write(line0)
        		counter_nop+=1
        #check for missing fields in called genotypes: (if find that length of format field is not hte same as genotype info field, and that hte genotype is not ./. then there's a problem)
        elif True in [(len(i.split(":"))!=len(myformat.split(":")) and i.split(":")[0]!="./.") for i in mygenoinfo]:
            errorVCF.write('# something is missing from genotype format for a called genotype \n')
            errorVCF.write(line0)
            counter_nop+=1
        # make sure that the way you calculated genotype counts checks out
        elif myCalled + myNoCalled != len(allCalls):
            errorVCF.write('# Something is wrong with genotype counts!\n')
            errorVCF.write(line0)
            counter_nop+=1
        #
        # check if there are more than 90% no-called genotypes -- made this very lenient just to get rid of egregious sites
        # really should get rid of them sooner.
        elif float(myNoCalled) > round(float(len(allCalls))* float(maxNoCallFrac)):
            errorVCF.write('# More than ' + str(maxNoCallFrac) + ' fraction of genotypes are  no-call at this site\n')
            errorVCF.write(line0)
            counter_nop+=1
        # check if AC score is correct; note that if there are no hets it will be a . not a 0
        # Note that AC is sum of 2*HomAlt + Het ; if there is no variation AC = . or is missing
        # remove sites where every individual is heterozygous (note you must do this AGAIN after you split by population)
        elif myHet==myCalled:
            errorVCF.write('# All genotypes are heterozygous\n')
            errorVCF.write(line0)
            counter_nop+=1
        # (no longer done here) AC, AN checks; not done because I know that GATK doesn't do AN/AC correctly, so I update them later in the script.
        #elif (myHet == 0 and myHomAlt==0) and 'AC' in infoFields and infoFields['AC']!=".":
        #    errorVCF.write('# Something is wrong with AC field-1!\n')
        #    errorVCF.write(line0)
        #    counter_nop+=1
        # if there are hets and/or HomAlt, and AC is in the infoField, then it should equal
        #elif (myHet != 0 or myHomAlt!=0) and infoFields['AC']!=myHet+2*myHomAlt:
        #    errorVCF.write('# Something is wrong with AC field-2!\n')
        #    errorVCF.write(line0)
        #    counter_nop+=1
        # check if AN is right:
        #elif (myCalled * 2) != infoFields['AN']:
        #    errorVCF.write('# Something is wrong with AN field!\n')
        #    errorVCF.write(line0)
        #    counter_nop+=1
        else:
            # want to update AN/AC for every line
            # lines that have no hets and no homs
            # note that I've made no changes to genotypes (unlike Clare's script)
            # if I had done that, then I'd need to recalculate myAN and myAC at the end here
            if myAC==0:
                # make a new line:
                mynewline=''
                infoFields['AC']="."
                infoFields['AF']="."
                infoFields['AN']=myCalled * 2
                mynewline += '\t'.join(line[:7]) # note: non-inclusive
                mynewline += '\t'
                mynewline += ';'.join('{0}={1}'.format(key, val) for key, val in sorted(infoFields.items()))
                mynewline += '\t'
                mynewline += '\t'.join(line[8:])
                mynewline += '\n'
                outVCF.write(mynewline)
                counter_p+=1
            # lines with variation:
            else:
                mynewline=''
                infoFields['AC']=myAC
                infoFields['AN']=myAN
                infoFields['AF']=round(float(myAC))/(float(myAN))
                mynewline += '\t'.join(line[:7])
                mynewline += '\t'
                mynewline += ';'.join('{0}={1}'.format(key, val) for key, val in sorted(infoFields.items()))
                mynewline += '\t'
                mynewline += '\t'.join(line[8:])
                mynewline += '\n'
                outVCF.write(mynewline)
                counter_p+=1
    # write out final information
    errorVCF.write('## This many lines failed a filter and werent printed' + '\t' + str(counter_nop) + '\n')
    errorVCF.write('## This many lines PASSED and were printed' + '\t' + str(counter_p) + '\n')
    errorVCF.close()
    outVCF.close()
    inVCF.close()


#### run the function ##########
main_vcf_check(filepath,outname,errorname)

# skipping a few things from clare's script
# not bothering to update PASS flag for no call genotypes
# not bothering to update alt allele for sites that become invariant to "." -- leaving as is
# but making sure the AC= .
sys.exit()
