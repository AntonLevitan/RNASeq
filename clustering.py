import pandas as pd
import matplotlib.pyplot as plt
from sklearn import cluster
import numpy as np

# data = pd.read_csv('drug_RNA_clustering.csv', index_col='Standard name')
drug_screen = pd.read_excel('FLC_network.xlsx', sheet_name='Drug Screen').drop('S score (FLC-NoDrug)', axis=1)

rna_seq = pd.read_excel('FLC_network.xlsx', sheet_name='RNA Seq')

df = drug_screen.merge(rna_seq).dropna().set_index('Standard name')

principalDf = df

k_means = cluster.KMeans(n_clusters=2, random_state=0)
final = k_means.fit(principalDf)

clusters = k_means.fit(principalDf).predict(principalDf)
centers = final.cluster_centers_

fig = plt.figure()
ax = fig.add_subplot(111)
scatter = ax.scatter(principalDf.iloc[:, 0], principalDf.iloc[:, 1], c=clusters, s=50)
for i, j in centers:
    ax.scatter(i, j, s=50, c='red', marker='+')
plt.xlabel(principalDf.columns.values[0])
plt.ylabel(principalDf.columns.values[1])
plt.colorbar(scatter)
plt.savefig('clustering_rna_drug2.png')
fig.show()


