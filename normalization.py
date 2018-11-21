import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def deseq(counts):
    # possible implementation of R library
    # need to validate
    # counts = counts[np.alltrue(counts, axis=1)]
    logcounts = np.log(counts)
    loggeommeans = np.mean(logcounts)
    sf = np.exp(np.median(logcounts - loggeommeans, axis=0))
    return sf


def zfkpm(fpkm):
    # https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-778
    # assume "fpkm" is a NumPy array of log2(fpkm) values.
    # looks like a good implementation
    # takes the right side of a gaussian dstribution
    # 'mirror' it to a full one
    # calculates z-score for each gene
    kernel = gaussian_kde(fpkm)
    xi = np.linspace(fpkm.min(), fpkm.max(), 100)
    yi = kernel.evaluate(xi)
    mu = xi[np.argmax(yi)]
    U = fpkm[fpkm > mu].mean()
    sigma = (U - mu) * np.sqrt(np.pi / 2)
    zFPKM = (fpkm - mu) / sigma
    # threshold of -3 should be used according to the paper
    return zFPKM


def rpk(counts, start, end):
    # normalizes read counts by gene length
    length = abs(np.array(start) - np.array(end))
    rpk = 1000 * np.array(length) / np.array(counts)
    return rpk


def rpm(counts):
    # normalizes by total library reads
    total = sum(counts)
    rpm = np.array(counts) / total
    return rpm


def rpkm(counts, start, end):
    # normalizes first by total and then by length
    rpm_counts = rpm(counts)
    rpkm = rpk(rpm_counts, start, end)
    return rpkm


def fpkm():
    # for pairwise alignment- fragments instead of reads, look it up when needed
    pass


def tpm(counts, start, end):
    # normalizes first by length and then by total
    rpk_counts = rpk(counts, start, end)
    tpm = rpm(rpk_counts)
    return tpm


data = pd.read_excel('A1_in_A_A17_in_A.xlsx')
x = list(data['ref_reads'])
logcounts = np.log(x)
print(logcounts)
print(zfkpm(logcounts))
