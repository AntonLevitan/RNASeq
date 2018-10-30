import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse

TRANSLATION_FILE = "translation table orf19--_A22 complete list.xlsx"
READS_COUNT = "ca25_count.txt"

moving_average_10_flag = '--ma10'
moving_average_100_flag = '--ma100'


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(moving_average_10_flag, dest='ma10', action="store_true")
    parser.add_argument(moving_average_100_flag, dest='ma100', action="store_true")

    return parser.parse_args()


def log2_reads(x):

        """Applying log2 to a number. Conditionally catching the error of log2(0)."""

        if x > 0:
            return math.log2(x)
        else:
            return 0


def plot_reads():

    """Reads a count of reads per feature .txt file.
    Loads a translation .xslx file (from orf19 to A22).
    Maps the reads on all the annotated features in A22.
    Transofrms the values into log2 scale.
    Takes a moving average of the log2 reads values of the features"""

    settings = arguments()
    features = pd.read_excel(TRANSLATION_FILE).drop("orf192", axis=1)
    reads = pd.read_csv(READS_COUNT, delimiter="\t", )
    reads.columns = ["orf19", "read_count"]
    mapped_reads = features.merge(reads, on="orf19")
    mapped_reads['chromosome'] = mapped_reads['A22'].str.split('_', 1).str[0]
    chromosomes = mapped_reads['chromosome'].unique()
    mapped_reads['log2_reads'] = mapped_reads['read_count'].apply(lambda x: log2_reads(x))
    mapped_reads['ma_log2_reads_10'] = mapped_reads['log2_reads'].rolling(window=10).mean()
    mapped_reads['ma_log2_reads_100'] = mapped_reads['log2_reads'].rolling(window=100).mean()

    for i in range(len(chromosomes)):
        mapped_reads['x'] = mapped_reads[mapped_reads['chromosome'] == chromosomes[i]]['log2_reads'].rolling(window=100).mean()    

    mapped_reads.to_excel("mapped_reads.xlsx")

    fig, axs = plt.subplots(1, 8, sharey=True, figsize=(20, 4))
    fig.subplots_adjust(hspace=.5, wspace=.001)
    axs = axs.ravel()

    for i in range(len(chromosomes)):
        chromosome = mapped_reads[mapped_reads['chromosome'] == chromosomes[i]]
        if settings.ma10:
            axs[i].bar(chromosome['A22'], chromosome["ma_log2_reads_10"], width=1.0)
        elif settings.ma100:
            axs[i].bar(chromosome['A22'], chromosome["ma_log2_reads_100"], width=1.0)
        else:
            axs[i].bar(chromosome['A22'], chromosome["log2_reads"], width=1.0)

    mapped_reads.to_excel("mapped_reads.xlsx")

    if settings.ma10:
        plt.savefig('ca25_ma_10_log2_reads.png')
    elif settings.ma100:
        plt.savefig('ca25_ma_100_log2_reads.png')
    else:
        plt.savefig('ca25_log2_reads.png')


if __name__ == "__main__":
    plot_reads()
