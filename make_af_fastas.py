# from pyrosetta import *
# init()
import os
import pandas as pd
df=pd.read_csv('baker_ting_energies_csts_sort.csv')
accept_list=list(df['name'])
# for files in os.listdir('interface_decoys_sorted/'):
# 	pose=pose_from_pdb('interface_decoys_sorted/'+files)
# 	with open('fastas/'+files.split('.')[0]+'.fa','w') as f:
# 		f.write('>'+files.split('.')[0]+'\n')
# 		f.write(pose.chain_sequence(1)+'\n')
names=[]
seqs=[]
with open('baker_ting_seqs.fa') as f:
	accept=0
	for x in f:
		if '>' in x:
			cur_name=x.replace('>','').strip()
			if cur_name in accept_list:
				accept=1
				names.append(cur_name)
		elif accept==1:
			seqs.append(x.strip())
			accept=0
for a, b in zip(names, seqs):
	with open('fastas_bt/'+a+'.fa','w') as f:
		f.write('>'+a+'\n')
		f.write(b+'\n')