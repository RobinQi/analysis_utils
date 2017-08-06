import os,sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

tmp = pd.read_table('pathway.txt',sep='\t')
fig, ax = plt.subplots()
y_pos = np.arange(len(tmp['KEGG']))
p_value= -np.log10(tmp['pValue'])
ax.barh(y_pos,p_value,color='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_yticks(y_pos+0.5)
ax.set_yticklabels(tmp['KEGG'])
ax.set_xlabel('-log10(pValue)')
ax.set_title('KEGG')

plt.show()