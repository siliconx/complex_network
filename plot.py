import os
import pandas as pd
import matplotlib.pyplot as plt

graphs = ['er', 'ws', 'ba']

for g in graphs:
    file = g + '.csv'
    if not os.path.exists(file):
        continue
    df = pd.read_csv(file, sep='\t')
    for v in ['w', 'q']:
        var = df[df['variable'] == v]
        new_var = var[['simula', 'formula']]
        new_var.index = var[v]
        new_var = new_var.sort_index()
        new_var.plot(title='%s_%s' % (g, v))
plt.show()
