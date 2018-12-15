import pandas as pd

location = pd.read_csv('/home/user/Desktop/RNASeq/Data/A22_current_features.gff', sep='\t', skiprows=25)
location = location[location.iloc[:, 2] == 'gene']
location['name'] = location.iloc[:, 8].str.split(';', 1).str[0].str.split('=', 1).str[1]
location = location[location['name'].str.split('_',  2).str[2] == 'A']
location['start'] = location.iloc[:, 3]
location['end'] = location.iloc[:, 4]
location['length'] = abs(location['start'] - location['end'])
print(location.head())
print(len(location['length']))
