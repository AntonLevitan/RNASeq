import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np

REFERENCE_READS = "CA_count.txt"
TRANSLATION_FILE = "translation table orf19--_A22 complete list.xlsx"
EXAMINED_READS = "ca25_count.txt"


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

    features = pd.read_excel(TRANSLATION_FILE).drop("orf192", axis=1)
    reads = pd.read_csv(read_count, delimiter="\t", )
    reads.columns = ["orf19", "read_count"]
    mapped_reads = features.merge(reads, on="orf19")
    mapped_reads['chromosome'] = mapped_reads['A22'].str.split('_', 1).str[0]
    mapped_reads['log2_reads'] = mapped_reads['read_count'].apply(lambda x: log2_reads(x))
    
    # TODO: 
    # calculate moving average for each chromosome independently
    mapped_reads['ma_log2_reads_10'] = mapped_reads['log2_reads'].rolling(window=10).mean()
    # reference['ma_log2_reads_100'] = reference['log2_reads'].rolling(window=100).mean()

    return mapped_reads


def matched_normalization():

    reference = mapping_reads(REFERENCE_READS)
    examined = mapping_reads(EXAMINED_READS)

    examined['ref_reads'] = reference['read_count']
    examined['ref_log2'] = reference['log2_reads']
    examined['ref_ma10'] = reference['ma_log2_reads_10']
    examined['reads_ratio'] = examined['ref_reads']/examined['read_count']
    examined['reads_log2_ratio'] = examined['ref_log2']/examined['log2_reads']
    examined['ma10_ratio'] = examined['ref_ma10']/examined['ma_log2_reads_10']
    examined['reads_ratio'] = examined['reads_ratio'].replace(np.inf, 0)
    examined['reads_log2_ratio'] = examined['reads_log2_ratio'].replace(np.inf, 0)
    examined['reads_ratio'] = examined['reads_ratio'].fillna(0)
    examined['reads_log2_ratio'] = examined['reads_log2_ratio'].fillna(0)
    examined.to_excel("mapped_reads.xlsx")

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
            axs[i].bar(chromosome['A22'], chromosome["ma10_ratio"], width=1.0)
        elif settings.ma100:
            axs[i].bar(chromosome['A22'], chromosome["ma_log2_reads_100"], width=1.0)
        else:
            axs[i].bar(chromosome['A22'], chromosome["reads_log2_ratio"], width=1.0)

    if settings.ma10:
        plt.savefig('ma10_ratio.png')
    elif settings.ma100:
        plt.savefig('ca25_ma_100_log2_reads.png')
    else:
        plt.savefig('reads_log2_ratio.png')


if __name__ == "__main__":
    plot_reads()
