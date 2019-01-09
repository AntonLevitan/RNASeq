import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import cluster

# classifier = pd.read_excel('FLC_network.xlsx', sheet_name='Classifier')
# classifier['Standard name'] = classifier['Standard name'].str[:-2]

# moursh = pd.read_excel('FLC_network.xlsx', sheet_name='Moursch')
# homann = pd.read_excel('FLC_network.xlsx', sheet_name='Homann')
# homann = homann.drop('Name', axis=1)
drug_screen = pd.read_excel('FLC_network.xlsx', sheet_name='Drug Screen').drop('S score (FLC-NoDrug)', axis=1)

rna_seq = pd.read_excel('FLC_network.xlsx', sheet_name='RNA Seq')
#
# df_list = [moursh, homann, drug_screen, rna_seq]
#
# df = df_list[0]
# for df_ in df_list[1:]:
#     df = df.merge(df_, on='Standard name')
#
# df.to_csv('network_FLC_combined.csv', index=False)

df = drug_screen.merge(rna_seq).dropna().set_index('Standard name')
# df.to_csv('drug2_RNA_clustering.csv')

# data = pd.read_csv('network_FLC_combined.csv', index_col='Standard name')
# norm_data = data
# norm_data['AVGRAD_moursch'] = np.log2(norm_data.AVGRAD_moursch)
# norm_data['aVGFoG_moursch'] = np.log2(norm_data.aVGFoG_moursch)
# norm_data['AVGRAD_homann'] = np.log2(norm_data.AVGRAD_homann)
# norm_data['aVGFoG_homann'] = np.log2(norm_data.aVGFoG_homann)
# norm_data = norm_data.round(2)
# x = norm_data.describe().round(2)
# norm_data.to_csv('log_network.csv')
# x.to_csv('log_network_stats.csv')

# df = pd.read_excel('alex_drug_screen.xlsx', index_col=0)

df['zscore_(FLC+Geld-FLC)'] = stats.zscore(list(df['S score (FLC+Geld-FLC)']))
df['pvalue (FLC+Geld-FLC)'] = stats.norm.sf(abs(df['zscore_(FLC+Geld-FLC)']))
df['zscore_RNA Seq A17 in vs. out'] = stats.zscore(list(df['RNA Seq A17 in vs. out']))
df['pvalue RNA Seq A17 in vs. out'] = stats.norm.sf(abs(df['zscore_RNA Seq A17 in vs. out']))
df2 = df[df['pvalue (FLC+Geld-FLC)'] < 0.05]
df2 = df2[df2['pvalue RNA Seq A17 in vs. out'] < 0.05]

principalDf = df2[df2['S score (FLC+Geld-FLC)'] < 0].iloc[:, :2]

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
