from pyrosetta import *
init()
import os
interface_indices=[]
from joey_utils import chain_selector, intergroup_selector, selector_intersection
for files in os.listdir('s219_decoys/'):
	pose=pose_from_pdb('s219_decoys/'+files)
	cha=chain_selector('A')
	chb=chain_selector('B')
	interface_selector=intergroup_selector(cha, chb)
	interface_cha=selector_intersection(cha, interface_selector)
	bools=interface_cha.apply(pose)
	ctr=8
	for a in bools:
		if a == 1 and ctr not in interface_indices:
			interface_indices.append(ctr)
		ctr+=1
print(interface_indices)
#combos:
#design 148, 170, 171, 177, 30, 218, 209, 211, 215 + interface
#fix 148, 170, 177, design interface
#fix 148, 170, 171, 177, 30, 218, 209, 211, 215, design interface

important_indices=[46, 81, 151, 219]
fix_partial=[46, 81, 151, 219, 171, 177, 176]
fix_full=[46, 81, 151, 219, 171, 177, 176, 211, 30, 209, 170, 218, 149, 214, 45, 148, 215]

full_design_str=''
partial_design_str=''
reduced_design_str=''

ctr1=0
ctr2=0
ctr3=0
for a in interface_indices:
	if a not in important_indices:
		full_design_str+=str(a-7)+' '
		ctr1+=1
	if a not in fix_partial:
		partial_design_str+=str(a-7)+' '
		ctr2+=1
	if a not in fix_full:
		reduced_design_str+=str(a-7)+' '
		ctr3+=1
print('full design')
print(full_design_str[:-1])
print('partial design')
print(partial_design_str[:-1])
print('reduced design')
print(reduced_design_str[:-1])
print(ctr1)
print(ctr2)
print(ctr3)