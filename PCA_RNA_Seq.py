import pandas as pd
import numpy as np
import random as rd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn import preprocessing

data = pd.read_excel('output.xlsx', header=1)
norm_data = data.loc[:'CR_10860C', 'A15_in_A:normalized.counts':'A3_out_C:normalized.counts'].dropna()
# print(norm_data.tail())
experiment_title = 'Alex Data'
k = 10

# generate example list genes
# genes = ['gene' + str(i) for i in range(1, 101)]


#  generate example treatments and samples
# wt = ['wt' + str(i) for i in range(1, 6)]
# ko = ['ko' + str(i) for i in range(1, 6)]


# data = pd.DataFrame(columns=[*wt, *ko], index=genes)
# # print(data.head())

# for gene in data.index:
#     data.loc[gene, 'wt1':'wt5'] = np.random.poisson(lam=rd.randrange(10, 1000), size=5)
#     data.loc[gene, 'ko1':'ko5'] = np.random.poisson(lam=rd.randrange(10, 1000), size=5)

# print(data.head())
# print(data.shape)

# similar to StandardScaler() + fit_transform
scaled_data = preprocessing.scale(norm_data.T)

# calling PCA object, fitting and transforming the scaled data to pca
pca = PCA()
pca.fit(scaled_data)
pca_data = pca.transform(scaled_data)

# % of variation of each pca
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]

# PCA per sample
pca_per_sample = pd.DataFrame(pca_data, index=[norm_data.columns.values], columns=labels)
print(pca_per_sample)

# take a look into the loading scores
# loading scores of PC1
loading_scores = pd.Series(pca.components_[0], index=norm_data.index.values)
# sort by absolute value
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
# take the top k genes; highest loading value=highest contribution to variance in a PC 
top_k_genes = sorted_loading_scores[0:k].index.values
# print(loading_scores[top_k_genes])


def scree_plot(variances, labels, title=experiment_title):

    ''' Charts the explained variance of each PC '''

    plt.bar(x=range(1, len(variances)+1), height=variances, tick_label=labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component')
    plt.title(title + ' Scree Plot')
    plt.savefig(title + ' Scree Plot' + '.png')  # not a typo in 'Scree'
    plt.show()


def pca_plot(pca_df, title=experiment_title):

    ''' Plots PC1 versus PC2 '''

    plt.scatter(pca_df.PC1, pca_df.PC2)

    for sample in pca_df.index:
        sample_label = str(sample[0]).split(':', 1)[0]
        plt.annotate(sample_label, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))

    plt.title(title + ' PCA Chart')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))
    plt.savefig(title + ' PCA Chart' + '.png')
    plt.show()


scree_plot(per_var, labels)
pca_plot(pca_per_sample)
