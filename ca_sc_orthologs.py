import pandas as pd

genes = pd.read_csv('top genes.csv')
orthologs = pd.read_csv('/home/user/Desktop/RNASeq/ca_sc_orthologs.csv')
orthologs['A22'] = orthologs['A22'].str[:-2]


top_orth = genes.merge(orthologs)
print(top_orth)
top_orth.to_csv('top100_sc_orth.csv')
