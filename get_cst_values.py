from pyrosetta import *
init('-run:preserve_header')
from pyrosetta.rosetta.protocols.constraint_generator import AddConstraints
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from joey_utils import apply_enzdes_constraints, tabulate_pdb_energy_lines
from pyrosetta.rosetta.core.scoring import ScoreType
import pandas as pd
import os
# features=['atom_pair_constraint','angle_constraint','dihedral_constraint']
# loop_indices=[]
# pose=pose_from_pdb('remarked_decoys/GKNE_msa_temp_1_fix_30_design_0.pdb')
# #pose=pose_from_pdb('remarked_decoys/QKTV_54_7_all_loop.pdb')
# #pose=pose_from_pdb('remarked_decoys/GKNE_loop_3_design_0.pdb')
# pose=apply_enzdes_constraints(pose, '5y4l_updated_4-25-23.cst')
# sf = get_fa_scorefxn()
# sf.set_weight(ScoreType.atom_pair_constraint, 1)
# sf.set_weight(ScoreType.coordinate_constraint, 1)
# sf.set_weight(ScoreType.angle_constraint, 1)
# sf.set_weight(ScoreType.dihedral_constraint, 1)
# print(sf(pose))
# #pose.dump_pdb('test.pdb')
# df=tabulate_pdb_energy_lines('test.pdb')
# df2=pd.DataFrame(df)
# df2.to_csv('test_msa.csv', index=False)
# for a in features:
# 	count=0
# 	print(a)
# 	for b, c in zip(df[a],df['pdb_number']):
# 		if b != 0 and type(c)==int:
# 			count+=1
# 			print(c)
#	print(count)
names=[]
atom_pair_constraint_sum=[]
angle_constraint_sum=[]
dihedral_constraint_sum=[]
coordinate_constraint_sum=[]
ctr=0
for files in os.listdir('meth_relax/'):
	# pose=pose_from_pdb('templates_relaxed/'+files)
	# pose=apply_enzdes_constraints(pose, '5y4l_updated_4-25-23.cst')
	# sf = get_fa_scorefxn()
	# sf.set_weight(ScoreType.atom_pair_constraint, 1)
	# sf.set_weight(ScoreType.coordinate_constraint, 1)
	# sf.set_weight(ScoreType.angle_constraint, 1)
	# sf.set_weight(ScoreType.dihedral_constraint, 1)
	# sf(pose)
	# pose.dump_pdb('templates_scored/'+files)
	df=tabulate_pdb_energy_lines('meth_relax/'+files)
	names.append(files.split('.')[0])
	atom_pair_constraint_sum.append(df['atom_pair_constraint'][1])
	angle_constraint_sum.append(df['angle_constraint'][1])
	dihedral_constraint_sum.append(df['dihedral_constraint'][1])
	coordinate_constraint_sum.append(df['coordinate_constraint'][1])
	ctr+=1
	print('$$$$$$$$$$$$')
	print(ctr)
	print('$$$$$$$$$$$$')
df2=pd.DataFrame(zip(names, atom_pair_constraint_sum, angle_constraint_sum, dihedral_constraint_sum, coordinate_constraint_sum), columns=['names','atom_pair_constraint_sum','angle_constraint_sum','dihedral_constraint_sum', 'coordinate_constraint_sum'])
df2.to_csv('meth_csts.csv',index=False)