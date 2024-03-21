from pyrosetta import *
init()
import os
all_seqs=[]
names=[]
for files in os.listdir('msa_sort/'):
	pose=pose_from_pdb('msa_sort/'+files)
	seq=pose.chain_sequence(1)
	seq_list=[]
	for a in seq:
		seq_list.append(a)
	print(seq_list[219-8])
	seq_list[219-8]='V'
	names.append(files.split('.')[0]+'_S219V')
	all_seqs.append(''.join(seq_list))
	seq_list[219-8]='N'
	names.append(files.split('.')[0]+'_S219N')
	all_seqs.append(''.join(seq_list))
with open('S219_muts.fa','w') as f:
	for a, b in zip(names, all_seqs):
		f.write('>'+a+'\n')
		f.write(b+'\n')
