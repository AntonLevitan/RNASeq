import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np

REFERENCE_READS = "CA_count_hapA.txt"
TRANSLATION_FILE = "translation table orf19--_A22 complete list.xlsx"
EXAMINED_READS = "A1_in_A_count"
PLOT_SUFFIX = '.png'
READ_COUNT_SUFFIX = '.txt'

READ_COUNT_FILE = EXAMINED_READS + READ_COUNT_SUFFIX

moving_average_10_flag = '--ma10'
moving_average_100_flag = '--ma100'


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(moving_average_10_flag, dest='ma10', action="store_true")
    parser.add_argument(moving_average_100_flag, dest='ma100', action="store_true")

    return parser.parse_args()


def log2_reads(x):

        """Applying log2 to a number. Conditionally catching the error of log2(0)"""

        if x > 0:
            return math.log2(x)
        else:
            return 0


def mapping_reads(read_count):

    # features = pd.read_excel(TRANSLATION_FILE).drop("orf192", axis=1)
    reads = pd.read_csv(read_count, delimiter="\t", header=None)
    reads = reads.iloc[:-5, :]
    reads.columns = ['A22', "read_count"]
    # reads['read_count'] = reads['read_count'].apply(lambda x: x + 1)
    # reads = features.merge(reads, on="orf19")
    reads['chromosome'] = reads['A22'].str.split('_', 1).str[0]
    reads['haplotype'] = reads['A22'].str.split('_', 2).str[2]
    reads = reads[reads['haplotype'] == 'A']
    reads['norm_reads'] = reads['read_count']/sum(reads['read_count'])*1000000
    print(sum(reads['read_count']))
    reads.to_excel('reads_hapA.xlsx')
    return reads


def matched_normalization():

    reference = mapping_reads(REFERENCE_READS)
    examined = mapping_reads(READ_COUNT_FILE)

    examined['ref_reads'] = reference['read_count']
    examined['norm_ref_reads'] = reference['norm_reads']
    examined['reads_ratio'] = examined['read_count']/examined['ref_reads']
    examined['norm_ratio'] = examined['norm_reads']/examined['norm_ref_reads']
    examined['reads_log2_ratio'] = examined['norm_ratio'].apply(lambda x: math.log2(x))
    #examined = examined.replace(0, 0.000001)
    #examined = examined.replace(np.inf, 0.000001)
    #examined = examined.fillna(0.000001)
    examined['log2_reads_ma10'] = examined['reads_log2_ratio'].rolling(window=10).mean()

    avg = examined.groupby(['chromosome']).mean()
    print(avg)
    examined.to_excel('reads.xlsx')

    return examined


def plot_reads():

    """Plots bar chart of reads in each chromosome"""

    settings = arguments()
    normalized = matched_normalization()
    chromosomes = normalized['chromosome'].unique()

    fig, axs = plt.subplots(1, 8, sharey=True, figsize=(20, 4))
    fig.subplots_adjust(hspace=.5, wspace=.001)
    axs = axs.ravel()

    for i in range(len(chromosomes)):
        chromosome = normalized[normalized['chromosome'] == chromosomes[i]]
        if settings.ma10:
            axs[i].bar(chromosome['A22'], chromosome["log2_reads_ma10"], width=1.0)
            axs[i].axhline(chromosome["log2_reads_ma10"].mean(), color='green', linewidth=2)
        elif settings.ma100:
            axs[i].bar(chromosome['A22'], chromosome["ma_log2_reads_100"], width=1.0)
        else:
            axs[i].bar(chromosome['A22'], chromosome["reads_log2_ratio"], width=1.0)
            axs[i].axhline(chromosome["log2_reads_ma10"].mean(), color='green', linewidth=2)

    if settings.ma10:
        plt.savefig(EXAMINED_READS + '_ma10_' + PLOT_SUFFIX)
    elif settings.ma100:
        plt.savefig(EXAMINED_READS + '_ma100_' + PLOT_SUFFIX)
    else:
        plt.savefig(EXAMINED_READS + '_' + PLOT_SUFFIX)


if __name__ == "__main__":
    mapping_reads(REFERENCE_READS)
