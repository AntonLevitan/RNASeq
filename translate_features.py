import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np

REFERENCE_READS = "CA_count.txt"
TRANSLATION_FILE = "translation table orf19--_A22 complete list.xlsx"
EXAMINED_READS = "ca25_count"
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
            return 0.000001


def mapping_reads(read_count):

    features = pd.read_excel(TRANSLATION_FILE).drop("orf192", axis=1)
    reads = pd.read_csv(read_count, delimiter="\t", )
    reads.columns = ["orf19", "read_count"]
    mapped_reads = features.merge(reads, on="orf19")
    mapped_reads['chromosome'] = mapped_reads['A22'].str.split('_', 1).str[0]

    return mapped_reads


def matched_normalization():

    reference = mapping_reads(REFERENCE_READS)
    examined = mapping_reads(READ_COUNT_FILE)

    examined['ref_reads'] = reference['read_count']
    examined['reads_ratio'] = examined['read_count']/examined['ref_reads']
    examined['reads_log2_ratio'] = examined['reads_ratio'].apply(lambda x: log2_reads(x))
    examined = examined.replace(0, 0.000001)
    examined = examined.replace(np.inf, 0.000001)
    examined = examined.fillna(0.000001)
    examined['log2_reads_ma10'] = examined['reads_log2_ratio'].rolling(window=10).mean()

    avg = examined.groupby(['chromosome']).mean()
    print(avg)
    examined.to_excel('mapped_reads.xlsx')

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

    if settings.ma10:
        plt.savefig(EXAMINED_READS + '_ma10_' + PLOT_SUFFIX)
    elif settings.ma100:
        plt.savefig(EXAMINED_READS + '_ma100_' + PLOT_SUFFIX)
    else:
        plt.savefig(EXAMINED_READS + '_' + PLOT_SUFFIX)


if __name__ == "__main__":
    plot_reads()
