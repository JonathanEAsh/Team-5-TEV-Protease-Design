import os
seqs=[]
names=[]
temps=[1,2,3]
fix=[30,50,70]
for i in temps:
	for j in fix:
		for files in os.listdir('mpnn/temp_'+str(i)+'/fix_'+str(j)+'/seqs/'):
			with open('mpnn/temp_'+str(i)+'/fix_'+str(j)+'/seqs/'+files) as f:
				line_ctr=1
				design_ctr=0
				for x in f:
					if line_ctr>2:
						if '>' not in x:
							seqs.append(x.strip())
							names.append('tev_msa_temp_'+str(i)+'_fix_'+str(j)+'_design_'+str(design_ctr))
							design_ctr+=1
					line_ctr+=1
seqs_nr=[]
names_nr=[]
for a, b in zip(names, seqs):
	if b not in seqs_nr:
		seqs_nr.append(b)
		names_nr.append(a)
print(len(seqs))
print(len(seqs_nr))
with open('tev_msa_mpnn_seqs.fa','w') as f:
	for a, b in zip(names_nr, seqs_nr):
		f.write('>'+a+'\n')
		f.write(b+'\n')