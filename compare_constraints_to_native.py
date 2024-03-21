import pandas as pd
names=[]
atom_diff=[]
angle_diff=[]
di_diff=[]
coord_diff=[]
df1=pd.read_csv('meth_csts.csv')
df2=pd.read_csv('tev_native_csts.csv')

#columns=df['names'],df['atom_pair_constraint_sum'],df['angle_constraint_sum'],df['dihedral_constraint_sum']

for a, b, c, d, e in zip(df1['names'],df1['atom_pair_constraint_sum'],df1['angle_constraint_sum'],df1['dihedral_constraint_sum'], df1['coordinate_constraint_sum']):
	motif=a.split('_')[0]
	names.append(a)
	for i, j, k, l, m in zip(df2['names'],df2['atom_pair_constraint_sum'],df2['angle_constraint_sum'],df2['dihedral_constraint_sum'], df2['coordinate_constraint_sum']):
		atom_diff.append(round(b-j,3))
		angle_diff.append(round(c-k,3))
		di_diff.append(round(d-l,3))
		coord_diff.append(round(e-m, 3))
df3=pd.DataFrame(zip(names, atom_diff, angle_diff, di_diff, coord_diff,), columns=['names', 'atom_diff', 'angle_diff', 'di_diff', 'coord_diff'])
df3.to_csv('meth_cst_diffs.csv', index=False)