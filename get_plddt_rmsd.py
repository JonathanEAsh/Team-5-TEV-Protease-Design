from pyrosetta import *
init()
import os
from joey_utils import get_b_factor, get_rmsd, chain_selector
import pandas as pd
rmsd=[]
plddt=[]
names=[]
for files in os.listdir('alt_out/best_models/'):
	pose=pose_from_pdb('alt_out/best_models/'+files)
	ros_pose=pose_from_pdb('baker_ting_decoys/'+files)
	cha=chain_selector('A')
	temp=0
	for i in range(len(pose.sequence())):
		temp+=get_b_factor(pose, i+1)
	plddt.append(round(temp/len(pose.sequence()),3))
	rmsd.append(round(get_rmsd(pose, ros_pose, cha, cha, cha, cha),3))
	names.append(files.split('.')[0])
df=pd.DataFrame(zip(names, plddt, rmsd), columns=['name','plddt','rmsd'])
df.to_csv('baker_ting_af_stats_sort_alt.csv', index=False)
