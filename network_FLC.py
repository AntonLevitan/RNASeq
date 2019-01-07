import pandas as pd

# classifier = pd.read_excel('FLC_network.xlsx', sheet_name='Classifier')
# classifier['Standard name'] = classifier['Standard name'].str[:-2]

moursh = pd.read_excel('FLC_network.xlsx', sheet_name='Moursch')
homann = pd.read_excel('FLC_network.xlsx', sheet_name='Homann')
homann = homann.drop('Name', axis=1)
drug_screen = pd.read_excel('FLC_network.xlsx', sheet_name='Drug Screen')
rna_seq = pd.read_excel('FLC_network.xlsx', sheet_name='RNA Seq')

df_list = [moursh, homann, drug_screen, rna_seq]

df = df_list[0]
for df_ in df_list[1:]:
    df = df.merge(df_, on='Standard name')

df.to_csv('network_FLC_combined.csv', index=False)
