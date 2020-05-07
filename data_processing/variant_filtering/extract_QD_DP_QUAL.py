import os
import sys
import argparse
import csv
import gzip

def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="This script will extract QUAL, DP and QD for all variant sites")

    parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. Path to the VCF file. Should be gzipped.")

    parser.add_argument(
            "--scaffold", required=True,
            help="REQUIRED. Full name of scaffold")

    parser.add_argument(
            "--outfile", required=True,
            help="REQUIRED. Path to the output file.")

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    outfile = open(args.outfile, "w")
    scaff=str(args.scaffold)
    outfile.write("QUAL\tDP\tQD\n")
    with gzip.open(args.VCF, "r") as VCF:
        for line in VCF:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                line = line.split("\t")
                scaffold= line[0]
                info_col = line[7]
                # if the scaffold isn't the scaff you've chosen, skip rest of loop
                if scaffold!=scaff:
                   continue
                # only work with lines that contain QD (snps)
                if "QD" in info_col:
                    myqual=line[5]
                    myinfo=line[7]
                    # split info fields:
                    # instead of iterating through each one, make a dict.
                    infoFields=dict(s.split('=') for s in myinfo.split(";"))
                    myQD=infoFields["QD"]
                    myDP=infoFields["DP"]
                    outfile.write("\t".join(str(x) for x in (myqual, myDP, myQD)))
                    outfile.write("\n")
                else:
                    continue
    VCF.close()
    print("closing vcf now")
    outfile.close()
    print("closing outfile now")
sys.exit(main())

