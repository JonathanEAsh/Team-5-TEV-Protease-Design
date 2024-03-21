import os
names=[]
seqs=[]
for fix_dirs in os.listdir('mpnn_interface/'):
	for temp_dirs in os.listdir('mpnn_interface/'+fix_dirs+'/'):
		for files in os.listdir('mpnn_interface/'+fix_dirs+'/'+temp_dirs+'/seqs/'):
			with open('mpnn_interface/'+fix_dirs+'/'+temp_dirs+'/seqs/'+files) as f:
				line_ctr=1
				design_ctr=0
				for x in f:
					if line_ctr>2:
						if '>' not in x:
							seqs.append(x.strip())
							names.append(files.split('.')[0]+'_interface_'+fix_dirs+'_temp_'+temp_dirs[-1]+'_design_'+str(design_ctr))
							design_ctr+=1
					line_ctr+=1
names_nr=[]
seqs_nr=[]
for a, b in zip(names, seqs):
	if b not in seqs_nr:
		names_nr.append(a)
		seqs_nr.append(b)
print(len(seqs))
print(len(seqs_nr))
with open('tev_interface_seqs.fa','w') as f:
	for a, b in zip(names_nr, seqs_nr):
		f.write('>'+a+'\n')
		f.write(b+'\n')