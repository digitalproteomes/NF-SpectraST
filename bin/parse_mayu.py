#!/usr/bin/env python

import pandas as pd
import argparse
import sys


def main():
    """Extracts iProphet probability score threshold to achieve 1%FDR at the protein level."""
    parser = argparse.ArgumentParser(description='Extracts p threshold from Mayu output file.')
    parser.add_argument('mayuFile', metavar='MAYUFILE', help='The mayu output file')
    parser.add_argument('--level', help='protFDR or pepFDR (the level to use for FDR thresholding)')
    args = parser.parse_args()
    mayu_file_name = args.mayuFile
    level = args.level
    
    df = pd.read_csv(mayu_file_name)
    probability = df[df[level] < 0.01].tail(1)['IP/PPs']
    probability.reset_index(drop=True, inplace=True)
    if len(probability):
        sys.stdout.write('{}'.format(probability[0]))


if __name__ == "__main__":
    main()
