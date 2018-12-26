import pandas as pd
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import preprocessing
from scipy import stats

experiment_title = 'Alex Drug Screen'

data = pd.read_excel('alex_drug_screen.xlsx', index_col=0)
data = data.drop(labels=['Essintial', 'S score (FLC-NoDrug)', 'S score (FLC+Geld-FLC)'], axis=1)
# norm_data = data.loc[:'CR_10860C', 'A15_in_A:normalized.counts':'A3_out_C:normalized.counts'].dropna()
norm_data = data

k = 100
#
# # generate example list genes
# # genes = ['gene' + str(i) for i in range(1, 101)]
#
# #  generate example treatments and samples
# # wt = ['wt' + str(i) for i in range(1, 6)]
# # ko = ['ko' + str(i) for i in range(1, 6)]
#
# # data = pd.DataFrame(columns=[*wt, *ko], index=genes)
# # # print(data.head())
#
# # for gene in data.index:
# #     data.loc[gene, 'wt1':'wt5'] = np.random.poisson(lam=rd.randrange(10, 1000), size=5)
# #     data.loc[gene, 'ko1':'ko5'] = np.random.poisson(lam=rd.randrange(10, 1000), size=5)
#
# # similar to StandardScaler() + fit_transform
scaled_data = preprocessing.scale(norm_data.T)
#
# # calling PCA object, fitting and transforming the scaled data to pca
pca = PCA(n_components=2)
pca.fit(scaled_data)
pca_data = pca.transform(scaled_data)
#
# # % of variation of each pca
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]
#
# # PCA per sample
pca_per_sample = pd.DataFrame(pca_data, index=[norm_data.columns.values], columns=labels)
#
# # take a look into the loading scores
# # loading scores of PC1
loading_scores = pd.Series(pca.components_[0], index=norm_data.index.values)
# sort by absolute value
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
# take the top k genes; highest loading value=highest contribution to variance in a PC
top_k_genes = sorted_loading_scores[0:k].index.values
print(top_k_genes)
pd.Series(top_k_genes).to_csv('top genes.csv')


def scree_plot(variances, pc_labels, title=experiment_title):

    """ Charts the explained variance of each PC """

    plt.bar(x=range(1, len(variances)+1), height=variances, tick_label=pc_labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component')
    plt.title(title + ' Scree Plot')
    plt.savefig(title + ' Scree Plot' + '.png')  # not a typo in 'Scree'
    plt.show()


def pca_plot(pca_df, title=experiment_title):

    """ Plots PC1 versus PC2 """

    plt.scatter(pca_df.PC1, pca_df.PC2)

    for sample in pca_df.index:
        sample_label = str(sample[0]).split(':', 1)[0]
        plt.annotate(sample_label, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))

    plt.title(title + ' PCA Chart')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))
    plt.savefig(title + ' PCA Chart' + '.png')
    plt.show()


def pca_scatter_all(df):

    # features = df.columns
    features = ['S value No-Drug', 'S value Flu', 'S value FLC-Geld']
    # Separating out the features
    x = df.loc[:, features].values
    # Separating out the target
    y = df.index.values
    # Standardizing the features
    x = StandardScaler().fit_transform(x)


    # loading_values_plot(pca_per_sample)
    # scree_plot(per_var, labels)
    # pca_plot(pca_per_sample)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    # per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ['principal component 1', 'principal component 2'])

    finalDf = pd.concat([principalDf, pd.Series(df.index.values)], axis = 1)

    tested_feature_1 = 'principal component 1'
    tested_feature_2 = 'principal component 2'
    x = finalDf[tested_feature_1]
    y = finalDf[tested_feature_2]

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))
    axScatter = plt.axes(rect_scatter)
    plt.xlabel(tested_feature_1, fontsize=18)
    plt.ylabel(tested_feature_2, fontsize=18)

    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, s=6)

    xbins = np.linspace(min(x), max(x), 30)
    ybins = np.linspace(min(y), max(y), 30)

    axHistx.hist(x, bins=xbins)
    axHisty.hist(y, bins=ybins, orientation='horizontal')

    plt.savefig('both_' + tested_feature_1 + ' vs. ' + tested_feature_2 + '.png')
    plt.show()


df = pd.read_excel('alex_drug_screen.xlsx', index_col=0)
df = df[df['Essintial'].isnull()]
df = df.drop(labels='Essintial', axis=1)
df['zscore_(FLC-NoDrug)'] = stats.zscore(list(df['S score (FLC-NoDrug)']))
df['zscore_(FLC+Geld-FLC)'] = stats.zscore(list(df['S score (FLC+Geld-FLC)']))
df['pvalue (FLC-NoDrug)'] = stats.norm.sf(abs(df['zscore_(FLC-NoDrug)']))
df['pvalue (FLC+Geld-FLC)'] = stats.norm.sf(abs(df['zscore_(FLC+Geld-FLC)']))
df2 = df[(df['pvalue (FLC-NoDrug)'] < 0.05) | (df['pvalue (FLC+Geld-FLC)'] < 0.05)]

# df.to_csv('zscore_Alex.csv')
# df2.to_csv('zscore_Alex_filtered.csv')
#
# pca_scatter_all(df2)
# #
plt.scatter(df2['S score (FLC-NoDrug)'], df2['S score (FLC+Geld-FLC)'], s=6)
plt.show()
#
# pca_plot(pca_per_sample)
