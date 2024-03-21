names=[]
seqs=[]
with open('hhblits_align.fasta') as f:
	check=0
	for x in f:
		if '>' in x:
			names.append(x.strip())
			if check==1:
				seqs.append(curr_seq)
			else:
				check=1
			curr_seq=''
		else:
			curr_seq+=x.strip()
seqs.append(curr_seq)
print(len(names))
print(len(seqs))
with open('tev_align_format.fasta','w') as f:
	for a, b in zip(names, seqs):
		f.write(a+'\n')
		f.write(b+'\n')