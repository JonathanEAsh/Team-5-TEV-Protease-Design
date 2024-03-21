def get_id(seq1, seq2):
	indices_to_check=[]
	ctr=0
	for a in seq1:
		if a != '-':
			indices_to_check.append(ctr)
		ctr+=1
	match_ctr=0
	for a in indices_to_check:
		if seq1[a] == seq2[a]:
			match_ctr+=1
	return match_ctr / len(indices_to_check)
def get_coverage(seq):
	tev_seq=seqs[0]
	indices_to_check=[]
	ctr=0
	for a in tev_seq:
		if a != '-':
			indices_to_check.append(ctr)
		ctr+=1
	match_ctr=0
	coverage_ctr=0
	for a in indices_to_check:
		if seq[a] == tev_seq[a]:
			match_ctr+=1
		if seq[a] != '-':
			coverage_ctr+=1
	return match_ctr / len(indices_to_check), coverage_ctr / len(indices_to_check)

names=[]
seqs=[]
with open('tev_reduced_2.fasta') as f:
	for x in f:
		if '>' in x:
			names.append(x.replace('>','').strip())
		else:
			seqs.append(x.strip())

final_names=[]
final_seqs=[]
for a, b in zip(names, seqs):
	if a == 'tev':
		final_names.append(a)
		final_seqs.append(b)
	else:
		id_list=[]
		for c in final_seqs:
			if c != b:
				id_list.append(get_id(b,c))
		max_id=max(id_list)
		query_id, query_coverage = get_coverage(b)
		print(max_id)
		print(query_id)
		print(query_coverage)
		if max_id <= 0.9 and query_id >= 0.3 and query_coverage >= 0.5:
			print('accept')
			final_names.append(a)
			final_seqs.append(b)
		else:
			print('fail')
with open('tev_filtered_2.fasta','w') as f:
	for a, b in zip(final_names, final_seqs):
		f.write('>'+a+'\n')
		f.write(b+'\n')
