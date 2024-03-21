indices_to_check=[]
seqs=[]
names=[]
aas=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
with open('tev_align_format.fasta') as f:
	line_ctr=1
	for x in f:
		if '>' not in x:
			seqs.append(x.strip())
			if line_ctr==2:
				index_ctr=0
				for a in x.strip():
					if a in aas:
						indices_to_check.append(index_ctr)
					index_ctr+=1
		else:
			names.append(x.strip())
		line_ctr+=1
seqs_red=[]
print(indices_to_check)
for a in seqs:
	temp_seq=''
	for b in indices_to_check:
		temp_seq+=a[b]
	seqs_red.append(temp_seq)
with open('tev_reduced_2.fasta','w') as f:
	for a, b in zip(names, seqs_red):
		f.write(a+'\n')
		f.write(b+'\n')