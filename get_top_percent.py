import matplotlib.pyplot as plt
import numpy as np
import math
from pyrosetta import *
init()
from joey_utils import chain_selector, selector_intersection, intergroup_selector
seqs=[]
with open('tev_filtered_2.fasta') as f:
	line_ctr=1
	for x in f:
		if line_ctr==2:
			native_seq=x.strip()
		else:
			if '>' not in x:
				seqs.append(x.strip())
		line_ctr+=1
data=[]
resnum=8
for i in range(len(native_seq)):
	if native_seq[i]!='-':	
		curr_res=native_seq[i]
		match=0
		mismatch=0
		for a in seqs:
			if a[i] == curr_res:
				match+=1
			else:
				mismatch+=1
		#print(curr_res)
		data.append([resnum, match/(match+mismatch)])
		resnum+=1
interface_indices=[]
pose=pose_from_pdb('remarked_templates/tev_4.pdb')
cha=chain_selector('A')
chb=chain_selector('B')
interface=intergroup_selector(cha, chb, nearby_atom=7)
interface_protease=selector_intersection(interface, cha)
bools=interface_protease.apply(pose)
ctr=8
for i in bools:
	if i == 1:
		interface_indices.append(ctr)
	ctr+=1
data.sort(key=lambda x:x[1], reverse=True)
print(data)
_30_count=math.ceil(0.3*len(data))
_50_count=math.ceil(0.5*len(data))
_70_count=math.ceil(0.7*len(data))
fix_30=[]
fix_50=[]
fix_70=[]
ctr=1
for i in data:
	if ctr <= _30_count:
		fix_30.append(i[0])
	if ctr <= _50_count:
		fix_50.append(i[0])
	if ctr <= _70_count:
		fix_70.append(i[0])
	ctr+=1
for a in interface_indices:
	if a not in fix_30:
		fix_30.append(a)
	if a not in fix_50:
		fix_50.append(a)
	if a not in fix_70:
		fix_70.append(a)
#include S219?
key_indices=[46, 81, 151, 148, 170, 177]
for a in key_indices:
	if a not in fix_30:
		fix_30.append(a)
	if a not in fix_50:
		fix_50.append(a)
	if a not in fix_70:
		fix_70.append(a)
fix_30_paper=[3, 7, 9, 10, 11, 12, 14, 19, 25, 34, 36, 38, 42, 44, 46, 47, 48, 51, 52, 53, \
55, 61, 62, 64, 68, 81, 88, 89, 90, 92, 94, 100, 101, 103, 110, 113, 116, 117, \
126, 127, 129, 139, 140, 142, 143, 144, 146, 149, 151, 152, 154, 156, 160, 161, \
163, 165, 167, 169, 177, 186, 190, 198, 202, 211, 212, 221]
fix_50_paper=[2, 3, 7, 8, 9, 10, 11, 12, 13, 14, 21, 23, 25, 26, 27, 31, 32, 34, 35, 36, 37, \
38, 41, 42, 43, 44, 46, 47, 48, 51, 52, 53, 55, 59, 61, 62, 64, 68, 70, 72, 76, \
81, 85, 88, 89, 90, 91, 92, 93, 94, 95, 98, 100, 101, 103, 107, 109, 112, 113, \
115, 116, 117, 119, 123, 125, 126, 127, 129, 133, 134, 135, 139, 140, 141, 142, \
143, 144, 146, 147, 148, 149, 150, 151, 152, 153, 154, 156, 157, 160, 161, 163, \
165, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 182, 183, \
186, 190, 198, 200, 202, 204, 205, 208, 209, 211, 212, 213, 214, 215, 216, 217, \
218, 219, 220, 221]
fix_70_paper=[1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 21, 22, 23, 25, 26, 27, 31, \
32, 33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 51, 52, 53, \
55, 57, 59, 61, 62, 63, 64, 66, 68, 69, 70, 71, 72, 73, 76, 79, 80, 81, 83, 84, \
85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 100, 101, 103, 107, \
108, 109, 111, 112, 113, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, \
127, 129, 131, 133, 134, 135, 137, 139, 140, 141, 142, 143, 144, 145, 146, 147, \
148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 160, 161, 163, 164, 165, \
166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 182, 183, \
186, 187, 189, 190, 194, 196, 198, 200, 202, 203, 204, 205, 206, 207, 208, 209, \
211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221]
fix_30_no=[]
fix_50_no=[]
fix_70_no=[]
for a in fix_30:
	if a not in fix_30_paper:
		fix_30_no.append(a)
for a in fix_50:
	if a not in fix_50_paper:
		fix_50_no.append(a)
for a in fix_70:
	if a not in fix_70_paper:
		fix_70_no.append(a)
print(fix_30_no)
print(fix_50_no)
print(fix_70_no)
# print(fix_30)
# print(fix_50)
# print(fix_70)
# print(len(fix_30))
# print(len(fix_50))
# print(len(fix_70))
fix_30_str=''
fix_50_str=''
fix_70_str=''
for a in fix_30:
	fix_30_str+=str(a-7)+' '
for a in fix_50:
	fix_50_str+=str(a-7)+' '
for a in fix_70:
	fix_70_str+=str(a-7)+' '
# print(interface_indices)
# print('30%')
# print(fix_30_str)
# print('50%')
# print(fix_50_str)
# print('70%')
# print(fix_70_str)