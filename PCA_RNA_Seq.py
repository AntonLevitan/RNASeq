import pandas as pd
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import preprocessing
from scipy import stats
from sklearn import cluster

experiment_title = 'Combined Data'

data = pd.read_csv('network_FLC_combined.csv', index_col='Name')
norm_data = data.drop(labels=['Standard name'], axis=1)
norm_data['AVGRAD_moursch'] = np.log2(norm_data.AVGRAD_moursch)
norm_data['aVGFoG_moursch'] = np.log2(norm_data.aVGFoG_moursch)
norm_data['AVGRAD_homann'] = np.log2(norm_data.AVGRAD_homann)
norm_data['aVGFoG_homann'] = np.log2(norm_data.aVGFoG_homann)

corr = pd.DataFrame(np.corrcoef(norm_data.T), columns=norm_data.columns.values, index=norm_data.columns.values)

k = 10

scaled_data = preprocessing.scale(norm_data.T)
pca = PCA()
pca.fit(scaled_data)
pca_data = pca.transform(scaled_data)
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]
pca_per_sample = pd.DataFrame(pca_data, index=[norm_data.columns.values], columns=labels)
loading_scores = pd.Series(pca.components_[0], index=norm_data.index.values)
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
top_k_genes = sorted_loading_scores[0:k].index.values


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

    features = norm_data.columns.values
    x = df.loc[:, features].values
    y = df.index.values
    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'],
                               index=y)

    tested_feature_1 = 'principal component 1'
    tested_feature_2 = 'principal component 2'
    x = principalDf[tested_feature_1]
    y = principalDf[tested_feature_2]

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

    plt.savefig(tested_feature_1 + ' vs. ' + tested_feature_2 + '.png')
    plt.show()


def k_means_clustering_pca(df):

    features = norm_data.columns.values
    x = df.loc[:, features].values
    y = df.index.values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    tested_feature_1 = 'principal component 1'
    tested_feature_2 = 'principal component 2'
    principalDf = pd.DataFrame(data=principalComponents, columns=[tested_feature_1, tested_feature_2], index=y)
    k_means = cluster.KMeans(n_clusters=3, random_state=4)
    final = k_means.fit(principalDf)

    clusters = k_means.fit(principalDf).predict(principalDf)
    centers = final.cluster_centers_

    fig = plt.figure()
    ax = fig.add_subplot(111)
    scatter = ax.scatter(principalDf.iloc[:, 0], principalDf.iloc[:, 1], c=clusters, s=50)
    for i, j in centers:
        ax.scatter(i, j, s=50, c='red', marker='+')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))
    plt.colorbar(scatter)
    plt.savefig('clustering.png')
    fig.show()

    genes = pd.Series(clusters, index=y)

    return genes[genes == 1]


print(k_means_clustering_pca(norm_data))
# scree_plot(per_var, labels)
# pca_scatter_all(norm_data)
# pca_plot(pca_per_sample)
