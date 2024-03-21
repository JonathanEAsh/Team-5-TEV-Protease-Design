from pyrosetta import *
init('-run:preserve_header')
import os
from joey_utils import index_selector
from joey_utils import intergroup_selector
from joey_utils import make_task_factory
from joey_utils import make_move_map
from joey_utils import fast_relax_mover
from joey_utils import extract_pose_chain
from joey_utils import delete_selection_mover
from joey_utils import chain_selector
from joey_utils import total_energy
from joey_utils import pack_mover
from joey_utils import get_rmsd
from joey_utils import delete_selection_mover
from joey_utils import selector_union
from joey_utils import pack_mover, apply_coord_constraints
import argparse
from pyrosetta.rosetta.protocols.constraint_generator import AddConstraints
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from joey_utils import apply_enzdes_constraints, tabulate_pdb_energy_lines
from pyrosetta.rosetta.core.scoring import ScoreType
def calc_ddg_fold(native_pose, mutant_pose, score_function):
	ddg_fold=total_energy(mutant_pose, score_function)-total_energy(native_pose, score_function)
	return ddg_fold

def calc_ddg_bind(native_complex_pose, mutant_complex_pose, score_function):
		
	native_pdz_pose=extract_pose_chain(native_complex_pose, 'A')
	native_peptide_pose=extract_pose_chain(native_complex_pose, 'B')
	native_complex_energy=total_energy(native_complex_pose, score_function)
	native_peptide_energy=total_energy(native_peptide_pose, score_function)
	native_pdz_energy=total_energy(native_pdz_pose, score_function)
	dg_bind_native=native_complex_energy-(native_peptide_energy+native_pdz_energy)

	mutant_pdz_pose=extract_pose_chain(mutant_complex_pose, 'A')
	mutant_peptide_pose=extract_pose_chain(mutant_complex_pose, 'B')
	mutant_complex_energy=total_energy(mutant_complex_pose, score_function)
	mutant_peptide_energy=total_energy(mutant_peptide_pose, score_function)
	mutant_pdz_energy=total_energy(mutant_pdz_pose, score_function)
	dg_bind_mutant=mutant_complex_energy-(mutant_pdz_energy+mutant_peptide_energy)

	ddg_bind=dg_bind_mutant-dg_bind_native
	ddg_fold=mutant_complex_energy - native_complex_energy
	return ddg_bind, dg_bind_mutant, mutant_complex_energy, ddg_fold

def calc_rmsd(native_pose, mutant_pose):
	rmsd=get_rmsd(native_pose, mutant_pose, r_type='all')
	return rmsd

def calc_dg_bind_clamp(mutant_complex_pose, score_function):
	mutant_pdz_pose=extract_pose_chain(mutant_complex_pose, 'A')
	mutant_peptide_pose=extract_pose_chain(mutant_complex_pose, 'B')
	mutant_complex_energy=total_energy(mutant_complex_pose, score_function)
	tf=make_task_factory()

	mutant_peptide_pose=pack_mover(mutant_peptide_pose, score_function, tf)
	mutant_peptide_energy=total_energy(mutant_peptide_pose, score_function)
	
	mutant_pdz_pose=pack_mover(mutant_pdz_pose, score_function, tf)
	mutant_pdz_energy=total_energy(mutant_pdz_pose, score_function)
	dg_bind_mutant=mutant_complex_energy-(mutant_pdz_energy+mutant_peptide_energy)
	return dg_bind_mutant, mutant_pdz_energy


parser = argparse.ArgumentParser()
#parser.add_argument("--name", help="Design philosophy", type=str)
parser.add_argument("--og", help="Original protein", type=str)
parser.add_argument("--design_seq", help="Designed protein sequence", type=str)
#parser.add_argument("--design_num", help="Number of design", type=str)

#parser.add_argument("--design_num", help="Number of design", type=int)

args = parser.parse_args()

#name=args.name
og=args.og
design_seq=args.design_seq
#design_num=args.design_num
#design_num=args.design_num
#bb_num=og.split('_')[-1]
#base_name=og[:-2]
#og_file=og.split('_')[0]+'_'+og.split('_')[1]+'.pdb'
#design_num=og.split('_')[-1]
# motif=og.split('_')[0]
# for files in os.listdir('remarked_templates/'):
# 	if motif in files:
# 		pose_complex=pose_from_pdb('remarked_templates/'+files)
# 		break
#pose_complex=pose_from_pdb('remarked_templates/tev_4.pdb')
#pose_complex=pose_from_pdb('all_fix_bbs/'+og_file)
#mut=og[-1]
#base=og.replace('_S219'+mut,'')
#pose_complex=pose_from_pdb('msa_sort/'+base+'.pdb')
#base=og.split('_')[0]+'_'+og.split('_')[1]+'_'+og.split('_')[2]+'_'+og.split('_')[3]+'_'+og.split('_')[4]+'_'+og.split('_')[5]+'_'+og.split('_')[6]+'_'+og.split('_')[7]+'_'+og.split('_')[8]+'_'+og.split('_')[9]+'_'+og.split('_')[10]+'_'+og.split('_')[11]+'_'+og.split('_')[12]+'_'+og.split('_')[13]+'_'+og.split('_')[14]
pose_complex=pose_from_pdb('final/'+og+'.pdb')
og_seq=pose_complex.chain_sequence(1)
clamp_length=len(og_seq)
#pose_complex=pose_from_pdb('start/DLG2.pdb')
pose_isolate=extract_pose_chain(pose_complex,'A')
#pose_complex=extract_pose_chain(pose_complex, 'C')
#ch_c=chain_selector('C')
#pose_complex=delete_selection_mover(pose_complex, ch_c)
#pose_complex.dump_pdb('complex.pdb')
sfxn=get_fa_scorefxn()
mutant_index_list=[]
#complex_mutant_index_list=[]
mutation_list=[]
native_index_list=[]
#complex_native_index_list=[]
char_ctr=0
for chars in design_seq:
	if chars=='/':
		break
	elif chars!=og_seq[char_ctr]:
		mutant_index_list.append(char_ctr+1)
		#complex_mutant_index_list.append(char_ctr+1+len(pose_complex.chain_sequence(1)))
		mutation_list.append(chars)
	else:
		native_index_list.append(char_ctr+1)
		#complex_native_index_list.append(char_ctr+1+len(pose_complex.chain_sequence(1)))
	char_ctr+=1

mutated_indices=index_selector(mutant_index_list)
native_indices=index_selector(native_index_list)
shell_indices=intergroup_selector(mutated_indices, native_indices)
shell_bools=shell_indices.apply(pose_isolate)
#print(shell_bools)
shell_index_list=[]
shell_ctr=1
for x in shell_bools:
	if x == 1:
		shell_index_list.append(shell_ctr)
	shell_ctr+=1
mut_dict=dict(zip(mutant_index_list, mutation_list))
print(mut_dict)
#complex_mutated_indices=index_selector(complex_mutant_index_list)
#complex_native_indices=index_selector(complex_native_index_list)
#complex_shell_indices=intergroup_selector(complex_mutated_indices, complex_native_indices)
#complex_shell_bools=complex_shell_indices.apply(pose_complex)
#print(complex_shell_bools)
#complex_shell_index_list=[]
#shell_ctr=1
#for x in complex_shell_bools:
#	if x == 1:
#		complex_shell_index_list.append(shell_ctr)
#	shell_ctr+=1
#complex_mut_dict=dict(zip(complex_mutant_index_list, mutation_list))

#pdz_selector=chain_selector('A')
#peptide_selector=chain_selector('B')
#interface_selector=intergroup_selector(pdz_selector, peptide_selector)

mv=make_move_map(bb=shell_index_list, chi=shell_index_list, jump=shell_index_list)
#mv_complex=make_move_map(bb=complex_shell_index_list, chi=complex_shell_index_list, jump=complex_shell_index_list)
#tf_mutant_isolate=make_task_factory(None, shell_indices, None, mut_dict, ex12=False)
#tf_mutant_complex=make_task_factory(None, complex_shell_indices, ace2_selector, complex_mut_dict, ex12=False)
#tf_native_isolate=make_task_factory(None, shell_indices, None, None, ex12=False)
#tf_native_complex=make_task_factory(None, complex_shell_indices, ace2_selector, None, ex12=False)
tf_mutant=make_task_factory(None, shell_indices, None, mut_dict, ex12=False)
tf_native=make_task_factory(None, shell_indices, None, None, ex12=False)
#tf_repack_interface=make_task_factory(None, interface_selector, None, None)
#tf_mutant_clamp=make_task_factory(None, shell_indices, None, mut_dict)
#ddg_bind_list=[]
dg_bind_list=[]
dg_fold_list=[]
avg_dg_fold_list=[]
#ddg_fold_list=[]
#rmsd_list=[]
mut_pose_complex=Pose()
mut_pose_isolate=Pose()
best_energy=10000
for i in range(5):
	mut_pose_complex.detached_copy(pose_complex)
	sfxn=get_fa_scorefxn()
	mut_pose_complex=apply_coord_constraints(mut_pose_complex)
	mut_pose_complex=apply_enzdes_constraints(mut_pose_complex, 'cstfile_oydv.cst')
	sfxn.set_weight(ScoreType.atom_pair_constraint, 1)
	sfxn.set_weight(ScoreType.coordinate_constraint, 1)
	sfxn.set_weight(ScoreType.angle_constraint, 1)
	sfxn.set_weight(ScoreType.dihedral_constraint, 1)
	#fr_mutant_isolate=fast_relax_mover(score_function=sfxn, task_factory=tf_mutant_isolate, movemap=mv)
	#fr_mutant_complex=fast_relax_mover(score_function=sfxn, task_factory=tf_mutant_complex, movemap=mv_complex)
	#fr_native_isolate=fast_relax_mover(score_function=sfxn, task_factory=tf_native_isolate, movemap=mv)
	#fr_native_complex=fast_relax_mover(score_function=sfxn, task_factory=tf_native_complex, movemap=mv_complex)
	fr_mutant=fast_relax_mover(score_function=sfxn, task_factory=tf_mutant)
	#fr_native=fast_relax_mover(score_function=sfxn, task_factory=tf_native, movemap=mv)
	fr=fast_relax_mover(score_function=sfxn)
	#fr_clamp=fast_relax_mover(score_function=sfxn, task_factory=tf_mutant_clamp)
	#mut_pose_isolate.detached_copy(pose_isolate)
	
	
	#fr_mutant_isolate.apply(mut_pose_isolate)
	#fr_mutant_complex.apply(mut_pose_complex)
	#fr_native_isolate.apply(pose_isolate)
	#fr_native_complex.apply(pose_complex)
	#fr_mutant.apply(mut_pose_isolate)
	#fr_mutant.apply(mut_pose_complex)
	#fr_native.apply(pose_isolate)
	#fr_native.apply(pose_complex)
	fr_mutant.apply(mut_pose_complex)

	#mut_pose_complex=pack_mover(pose=mut_pose_complex, score_function=sfxn, task_factory=tf_repack_interface)
	#pose_complex=pack_mover(pose=pose_complex, score_function=sfxn, task_factory=tf_repack_interface)

	#mut_pose_isolate.dump_pdb('mut_pose_isolate.pdb')
	#pose_isolate.dump_pdb('nat_pose_isolate.pdb')
	#pose_complex.dump_pdb('nat_pose_complex.pdb')

	fr.apply(mut_pose_complex)
	#ddg_bind, dg_bind, dg_fold, ddg_fold=calc_ddg_bind(pose_complex, mut_pose_complex, sfxn)
	dg_bind, dg_fold = calc_dg_bind_clamp(mut_pose_complex, sfxn)
	#ddg_bind=calc_ddg_bind(pose_complex, mut_pose_complex, sfxn)
	#rmsd=calc_rmsd(pose_isolate, mut_pose_isolate)
	#rmsd_list.append(rmsd)
	#ddg_bind_list.append(ddg_bind)
	dg_bind_list.append(dg_bind)
	dg_fold_list.append(dg_fold)
	avg_dg_fold_list.append(dg_fold/clamp_length)
	#ddg_fold_list.append(ddg_fold)
	tot=total_energy(mut_pose_complex, sfxn)
	if tot < best_energy:
		best_energy=tot
		mut_pose_complex.dump_pdb('meth_relax/'+og+'.pdb')


#ace2_native_seq=pose_complex.chain_sequence(1)
#ace2_mutant_seq=mut_pose_complex.chain_sequence(1)
#ctr=1
#print('ACE2')
#for a, b in zip(ace2_native_seq, ace2_mutant_seq):
#	if a!=b:
#		print(a+str(ctr)+b)
#	ctr+=1
#rbd_native_seq=pose_complex.chain_sequence(2)
#rbd_mutant_seq=mut_pose_complex.chain_sequence(2)
#ctr=1
#print('RBD')
#for a, b in zip(rbd_native_seq, rbd_mutant_seq):
#	if a!=b:
#		print(a+str(ctr)+b)
#	ctr+=1
#
#rbd_native_seq=pose_isolate.sequence()
#rbd_mutant_seq=mut_pose_isolate.sequence()
#ctr=1
#print('isolate')
#for a, b in zip(rbd_native_seq, rbd_mutant_seq):
#	if a!=b:
#		print(a+str(ctr)+b)
#	ctr+=1
true_dg_fold=min(dg_fold_list)
true_dg_bind=dg_bind_list[dg_fold_list.index(true_dg_fold)]
true_avg_dg_fold=avg_dg_fold_list[dg_fold_list.index(true_dg_fold)]
#true_ddg_fold=min(ddg_fold_list)
#true_dg_fold=dg_fold_list[dg_bind_list.index(true_dg_bind)]
#true_dg_bind=dg_bind_list[ddg_fold_list.index(true_ddg_fold)]
#true_ddg_bind=ddg_bind_list[ddg_fold_list.index(true_ddg_fold)]
#true_rmsd=min(rmsd_list)
#name=og+'_design_'+str(design_num)
with open('relax_output/'+og+'.txt','w') as f:
	f.write(str(true_dg_bind)+','+str(true_dg_fold)+','+str(true_avg_dg_fold))