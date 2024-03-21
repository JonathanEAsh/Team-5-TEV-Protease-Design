seqs=[]
names=[]
es=[4,10,30,50]
for vals in es:
	with open('hhblits_e'+str(vals)+'.a3m') as f:
		accept=0
		for x in f:
			if '>' in x:
				if x.replace('>','').strip() not in names:
					names.append(x.replace('>','').strip())
					accept=1
			elif accept==1:
				accept=0
				seqs.append(x.strip().replace('-',''))
print(len(seqs))
print(seqs[0:5])
with open('hhblits_format.fasta','w') as f:
	for a, b in zip(names, seqs):
		f.write('>'+a+'\n')
		f.write(b+'\n')
