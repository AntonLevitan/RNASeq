import matplotlib.pyplot as plt
import argparse
import pandas as pd

moving_average_10_flag = '--ma10'
moving_average_100_flag = '--ma100'

PLOT_SUFFIX = '.png'
DATA_SUFFIX = '.xlsx'


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(moving_average_10_flag, dest='ma10', action="store_true")
    parser.add_argument(moving_average_100_flag, dest='ma100', action="store_true")

    return parser.parse_args()


def plot_reads():

    """Plots bar chart of reads in each chromosome"""
    
    normalized = pd.read_excel('7-32-A-1_S6_L001_R1_001_198-32-A-1_S28_L003_R1_001.xlsx')

    settings = arguments()
    chromosomes = normalized['chromosome'].unique()

    fig, axs = plt.subplots(1, len(chromosomes), sharey=True, figsize=(20, 4))
    fig.subplots_adjust(hspace=.5, wspace=.001)
    axs = axs.ravel()

    for i in range(len(chromosomes)):
        chromosome = normalized[normalized['chromosome'] == chromosomes[i]]
        chromosome['ma10'] = chromosome['tpm_ratio'].rolling(window=10).mean()
        chromosome['ma100'] = chromosome['tpm_ratio'].rolling(window=100).mean()
        if settings.ma10:
            axs[i].bar(chromosome['id'], chromosome['ma10'], width=1.0)
            axs[i].axhline(chromosome['ma10'].mean(), color='green', linewidth=2)
            plt.ylim(0.3, 1.8)
        elif settings.ma100:
            axs[i].bar(chromosome['id'], chromosome['ma100'], width=1.0)
            axs[i].axhline(chromosome['ma100'].mean(), color='green', alpha=0.7, linewidth=2)
            axs[i].axhline(0.5, color='red', linewidth=2, alpha=0.3)
            # axs[0].text('1 copy', 0.5, "{:.0f}".format(0.5), color="red", ha="right", va="center")    
            axs[i].axhline(1, color='red', linewidth=2, alpha=0.3)
            # axs[0].text('2 copies', 1, "{:.0f}".format(0.5), color="red", ha="right", va="center")
            axs[i].axhline(1.5, color='red', linewidth=2, alpha=0.3)
            # axs[0].text('3 copies', 1.5, "{:.0f}".format(0.5), color="red", ha="right", va="center")
            # plt.ylim(0.75, 1.7)
            axs[i].set_xlabel(chromosomes[i])
        else:
            axs[i].bar(chromosome['id'], chromosome['tpm_ratio'], width=1.0)
            axs[i].axhline(chromosome['tpm_ratio'].mean(), color='green', linewidth=2)

    if settings.ma10:
        plt.savefig('7-32-A-1_S6_L001_R1_001_198-32-A-1_S28_L003_R1_001' + PLOT_SUFFIX)
    elif settings.ma100:
        plt.savefig('7-32-A-1_S6_L001_R1_001_198-32-A-1_S28_L003_R1_001' + PLOT_SUFFIX)
    else:
        plt.savefig('7-32-A-1_S6_L001_R1_001_198-32-A-1_S28_L003_R1_001' + PLOT_SUFFIX)


plot_reads()
