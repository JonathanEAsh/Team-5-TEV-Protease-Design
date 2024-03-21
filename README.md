This repository contains all python scripts used to redesign the TEV protease. 
------------------------------------------------------------------------------
add_baker_ting.py - Compares interface sequences to full sequences provided by the Baker and Ting labs, installing mutations outside of the interface

all_tev_design_stats.csv - Design statistics of the 94 final design candidates, including interaction energy differences, constraint score differences, plddt, and rmsd

compare_constraints_to_native.py - Subtract native contraint values from design constraint values

compile_mpnn_seqs_interface.py - Compile MPNN-generated sequences from the interface round into one fasta file

compile_mpnn_seqs_msa.py - Compile MPNN-generated sequences from the MSA round into one fasta file

cstfile_oydv.cst - Constraints for the TEV protease

filter_fasta.py - Filter MSA based on query coverage, percent ID to query, and overall percent ID redundancy

format_fasta.py, format_for_logo.py - Format alignment by concatenating sequences to one line and removing query gaps

get_cst_values.py - Calculate constraint score sums for designs

get_final_alignment.py - Combine HHBlits alignments at different e-values, checking for redundancy

get_interface_indices.py - Return interface indices for MPNN to design based on PGCN node importance

get_plddt_rmsd.py - Compute pLDDT and RMSD for AlphaFold-modeled designs

get_top_percent.py - Return most conserved regions of the protease based on MSA percent ids for thermostabilization with MPNN

make_af_fastas.py - Read through final design set to write fastas for subsequent AlphaFold prediction

make_s219_muts.py - Install S219V/N point mutations on filtered MSA designs

mpnn_wrapper.py - Python wrapper for running ProteinMPNN

score_designs.py, score_interface.py, relax_submit.py - Reads through fasta file, grafting each sequence onto a corresponding structure in parallel with PyRosetta FastRelax

score_new_msa.py - Compute interaction energy differences to the native along the peptide for all designs

sort.py - Filter down designs based on constraints and interaction energies
