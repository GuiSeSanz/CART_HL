import pandas as pd
import pickle
import os

data_folder = os.path.join(os.getcwd(), 'Data')
files = list(filter(lambda x: x.endswith('_AUCs_filtered_BIC.pickle') , os.listdir(data_folder)))


files = ['CART_HighLow_L10.01_L20.1_AUCs_filtered_BIC.pickle']

for fl in files:
    a = pd.read_pickle(os.path.join(data_folder, fl))
    for k in a.keys():
        new_file = data_folder + '/' + fl.replace('.pickle', '_'+ str(k) +'_BIS.csv')
        print(new_file)
        a[k] = a[k].fillna(0)
        a[k].to_csv(new_file, sep='\t', header = True, index=True )



