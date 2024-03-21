import os
from pyrosetta import *
init()
from joey_utils import calc_single_res_intEs, chain_selector
import pandas as pd
names=[]
p7=[]
p6=[]
p5=[]
p4=[]
p3=[]
p2=[]
p1=[]
p1_=[]
p2_=[]
percent_muts=[]
og_pose=pose_from_pdb('template_relaxed/tev.pdb')
og_seq=og_pose.chain_sequence(1)
cha=chain_selector('A')
chb=chain_selector('B')
sfxn=get_fa_scorefxn()
og_energies=calc_single_res_intEs(og_pose, sfxn, chb, cha)
print(og_energies)
for files in os.listdir('meth_relax/'):
	pose=pose_from_pdb('meth_relax/'+files)
	energies=calc_single_res_intEs(pose, sfxn, chb, cha)
	print(energies)
	print(len(energies))
	seq=pose.chain_sequence(1)
	ctr=0
	for a, b in zip(seq, og_seq):
		if a != b:
			ctr+=1
	names.append(files.split('.')[0])
	p7.append(round(energies[0]-og_energies[0],3))
	p6.append(round(energies[1]-og_energies[1],3))
	p5.append(round(energies[2]-og_energies[2],3))
	p4.append(round(energies[3]-og_energies[3],3))
	p3.append(round(energies[4]-og_energies[4],3))
	p2.append(round(energies[5]-og_energies[5],3))
	p1.append(round(energies[6]-og_energies[6],3))
	p1_.append(round(energies[7]-og_energies[7],3))
	p2_.append(round(energies[8]-og_energies[8],3))
	percent_muts.append(round(ctr/len(og_seq),3))
df=pd.DataFrame(zip(names, p7, p6, p5, p4, p3, p2, p1, p1_, p2_, percent_muts), columns=['name','p7','p6','p5','p4','p3','p2','p1','p1*','p2*','percent_muts'])
df.to_csv('meth_energies.csv',index=False)