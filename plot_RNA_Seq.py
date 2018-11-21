import matplotlib.pyplot as plt
import argparse
import pandas as pd

moving_average_10_flag = '--ma10'
moving_average_100_flag = '--ma100'

PLOT_SUFFIX = '.png'
DATA_SUFFIX = '.xlsx'

# 'A1_in_A_A17_in_A' = input('Input examined filename here: ')


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(moving_average_10_flag, dest='ma10', action="store_true")
    parser.add_argument(moving_average_100_flag, dest='ma100', action="store_true")

    return parser.parse_args()


def plot_reads():

    """Plots bar chart of reads in each chromosome"""
    normalized = pd.read_excel('A1_in_A_A17_in_A.xlsx')

    settings = arguments()
    chromosomes = normalized['chromosome'].unique()

    fig, axs = plt.subplots(1, 8, sharey=True, figsize=(20, 4))
    fig.subplots_adjust(hspace=.5, wspace=.001)
    axs = axs.ravel()

    for i in range(len(chromosomes)):
        chromosome = normalized[normalized['chromosome'] == chromosomes[i]]
        if settings.ma10:
            axs[i].bar(chromosome['A22'], chromosome['log2_reads_ma10'], width=1.0)
            axs[i].axhline(chromosome["log2_reads_ma10"].mean(), color='green', linewidth=2)
        elif settings.ma100:
            axs[i].bar(chromosome['A22'], chromosome['log2_reads_ma100'], width=1.0)
        else:
            axs[i].bar(chromosome['A22'], chromosome['reads_log2_ratio'], width=1.0)
            axs[i].axhline(chromosome['reads_log2_ratio'].mean(), color='green', linewidth=2)

    if settings.ma10:
        plt.savefig('A1_in_A_A17_in_A' + '_ma10' + PLOT_SUFFIX)
    elif settings.ma100:
        plt.savefig('A1_in_A_A17_in_A' + '_ma100' + PLOT_SUFFIX)
    else:
        plt.savefig('A1_in_A_A1_in_B' + PLOT_SUFFIX)


plot_reads()
