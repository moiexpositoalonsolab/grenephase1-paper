import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
import math


parser = argparse.ArgumentParser()
parser.add_argument('-i',type = str, action= 'store',dest='input',help='the selection atlas file')
parser.add_argument('-s',type = str, action= 'store',dest='size',help='Chromosome size')
parser.add_argument('-o',type = str, action = 'store', dest = 'output',help = "the prefix of the output files")

args = parser.parse_args()

#read dataset


pvalues = pd.read_csv(args.input, sep = '\t',header = None)

r,c = pvalues.shape


ch_size = [100000]
space = 0
with open(args.size,"r") as f:
	for line in f:
		line = line.strip("\n")
		ch_size.append(int(line)+space)

ch_size_cum = np.cumsum(ch_size)

fig = plt.figure(figsize=[30,10],dpi=600)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(True)

colors = ['dimgrey','silver']
x_labels = []
x_labels_pos = []


for i in range(1,len(ch_size_cum)):
	x_labels_pos.append((ch_size_cum[i]-ch_size_cum[i-1])/2+ch_size_cum[i-1])


for i in range(len(ch_size) - 1 ):
	index = np.where(pvalues.values[:,0] == i+1)[0]
	locus_start = np.array(pvalues.values[index,1]) + ch_size_cum[i]
	locus_end = np.array(pvalues.values[index,2]) + ch_size_cum[i]
	x = (locus_start + locus_end ) / 2.0
	width = (locus_end - locus_start)/ 50000
	p = np.array(pvalues.values[index,3])
	col = colors[i % 2]
	plt.vlines(x,ymin=0,ymax=p,color = col,linewidth=width)

plt.axhline(y=5, color='r', linestyle='--')

ax.set_xticks(x_labels_pos)
ax.set_xticklabels(np.arange(1,len(ch_size_cum)))
ax.set_xlim([0, ch_size_cum[-1]])
ax.set_ylim([0, np.max(pvalues.values[:,3]) + 1 ])
plt.xlabel('Chromosome',fontsize=30)
plt.ylabel('Parallele Shift Occurrence',fontsize=30)
plt.xticks(fontsize = 25,fontweight = 600)
plt.yticks(fontsize = 25,fontweight = 600)

plt.savefig(args.output+".pdf",format = "pdf",orientation = "landscape")