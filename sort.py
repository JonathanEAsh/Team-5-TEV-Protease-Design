import pandas as pd
import subprocess
features_1=['p7','p6','p5','p4','p3','p2','p1','p1*','p2*','percent_muts']
features_1_compare=['p6','p5','p2']
features_2=['atom_diff','angle_diff','di_diff','coord_diff']
neg_feats=[]
df1=pd.read_csv('baker_ting_energies.csv')
df2=pd.read_csv('baker_ting_cst_diffs.csv')
#df2_og=pd.read_csv('../loop_msa_design_constraints_differences_redo.csv')
#df1=df1.sort_values(by='p2',axis=0,ignore_index=True)
ctr=0
accept=[]
for a in df1['name']:
	bad=0
	for b in features_1_compare:
		# if b == 'p6' or b == 'p5':
		# 	if df1[b][ctr]>1:
		# 		bad+=1
		if df1[b][ctr]>0:
			bad+=1
	ctr+=1
	ctr2=0
	#if bad <= 1:
	for b in df2['names']:
		if b == a:
			for c in features_2:
				if c=='coord_diff': 
					if df2[c][ctr2]>2:
						bad+=1
				elif df2[c][ctr2]>1:
					bad+=1
		ctr2+=1
	if bad<=1:
		accept.append(a)
	# else:
	# 	print(coord)
print(len(accept))
ctr=0
accept_ctr=0
df3=pd.DataFrame(columns=['name']+features_1+features_2)
for a in accept:
	temp=[a]
	ctr=0
	for b in df1['name']:
		if b == a:
			for c in features_1:
				temp.append(df1[c][ctr])
		ctr+=1
	ctr=0
	for b in df2['names']:
		if b == a:
			for c in features_2:
				temp.append(df2[c][ctr])
		ctr+=1
	# print(temp)
	df3=df3.append(pd.Series(temp, index=['name']+features_1+features_2), ignore_index=True)
# 	subprocess.call(['cp','interface_decoys/'+a+'.pdb','interface_decoys_sorted/'])
# # muts=[]
# # d4=pd.read_csv('new_msa_design_muts_temp.csv')
# # for a in d3['name']:
# # 	for b, c in zip(d4['name'],d4['num_muts']):
# # 		if b == a:
# # 			muts.append
df3.to_csv('baker_ting_energies_csts_sort_alt.csv',index=False)