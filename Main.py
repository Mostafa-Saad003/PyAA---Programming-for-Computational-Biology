# # This is a main script for testing all the modules!
import fetch as f, characterization as ch, categorization as ca, searching as s, structure_prediction as sp
import time
#
# # This is the data that we'll test against all the modules
# vdac = f.fetch_proteins("https://rest.uniprot.org/uniprotkb/P21796.fasta", "VDAC", 1)
# transmembrane_proteins = f.fetch_proteins("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28DNA-binding%20proteins%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28existence%3A1%29%20AND%20%28proteins_with%3A56%29%20AND%20%28annotation_score%3A5%29%20AND%20%28length%3A%5B1%20TO%20200%5D%29", "Transmembrane", 62)
# dna_binding = f.fetch_proteins("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29%20AND%20%28model_organism%3A9606%29%20AND%20%28proteins_with%3A1%29%20AND%20%28proteins_with%3A23%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28annotation_score%3A5%29%20AND%20%28length%3A%5B201%20TO%20400%5D%29", "DNA-Binding Proteins", 58)
# binding_protein_1 = f.fetch_proteins('https://rest.uniprot.org/uniprotkb/Q13541.fasta', '4E-binding protein I', 1)
# binding_protein_2 = f.fetch_proteins('https://rest.uniprot.org/uniprotkb/P70445.fasta', '4E-binding protein II', 1)
# dessication_protein = f.fetch_proteins('https://rest.uniprot.org/uniprotkb/Q8LAU8.fasta', 'Dessication-related protein', 1)
#
# # First, we test protein characterization. This script will generate 2 .txt file that contain the majority of functions
# # One file for a list of proteins, and another for a single protein
#
# ch.export(transmembrane_proteins, "Transmembrane_proteins", "transmembrane.txt", 62)
# ch.export(dna_binding, "DNA-binding_proteins", "dna_binding.txt", 58)
# ch.export(vdac, "VDAC", "VDAC.txt", 1)
# # ---------------------------------------------------------------------------------------------------------------------
#
# # Second, we test protein categorization according to one-of-two-groups
# # Note this method will take a long time, somewhere around 3-5 minutes, to run
# print("Categorization according to one-of-two-groups, by smallest distance:")
# print(ca.categorize_by_distance(vdac, transmembrane_proteins, dna_binding))
# print("-----------------------------------------------")
#
# # To test the accuracy of this method, take the first 10 proteins of the transmembrane protein list and test them against the same list except the first 10 proteins
# # If the accuracy is 100%, all of them should belong to transmembrane proteins
#
# testing_proteins = transmembrane_proteins[:10]
# counter = 0
# for protein in testing_proteins:
#     result = ca.categorize_by_distance(protein, transmembrane_proteins[11:], dna_binding)
#     if result == "This protein belongs to category 1":
#         counter += 1
#     else:
#         continue
#
# print("The accuracy is: ",(counter/10) * 100)
# print("-----------------------------------------------")
#
# # Testing categorization (natively folded vs. natively unfolded) according to mean hydrophobicity
# print('Natively Folded/Unfolded according to mean hydrophobicity:')
# print("4E-binding protein I")
# ca.folded_or_unfolded(binding_protein_1)
# print('\n4E-binding protein II')
# ca.folded_or_unfolded(binding_protein_2)
# print('\nDessication-related protein')
# ca.folded_or_unfolded(dessication_protein)
# print("-----------------------------------------------")
#
# # ---------------------------------------------------------------------------------------------------------------------
#
# # Third, we test secondary structure prediction
#
# # First, by Chou and Fasman
# print("Secondary structure prediction by Chou and Fasman's statistical propensity method:")
# print(sp.get_secondary_structures_by_chou_and_fasman(vdac))
# print("-----------------------------------------------")
#
# # Second, by dihedral angles
# coordinates = sp.getCoordinates('2jk4.pdb') # This is the pdb file of VDAC1
# phis_psis = sp.calculate_phi_psi(vdac, coordinates)
# print("Dihedral angles: ")
# print(phis_psis)
# print("-----------------------------------------------")
# sp.plot_ramachandran(phis_psis)
#
# # Third, by hydrophobicity profile plotting
# sp.plot_hydrophobicity_profile(vdac, 5)
#
# # ---------------------------------------------------------------------------------------------------------------------
#
# # Fourth and finally, we test pattern search algorithms
#
# protein_seq = "MAVPPTYYADLPPTYYGKSPPTA"
# pat = "PPT"
#
# print("Pattern searching using Naive algorithm:")
# start_time = time.time() * 1000000000
# print(s.Naive(pat, protein_seq))
# end_time = time.time() * 1000000000
# print("Naive takes ", start_time - end_time, " to run")
#
# print("Pattern searching using Boyer-Moore algorithm:")
# start_time2 = time.time() * 1000000000
# print(s.BoyerMoore(pat, protein_seq))
# end_time2 = time.time() * 1000000000
# print("Boyer Moore takes ", start_time2 - end_time2, " to run")
#
# print("Pattern searching using KMP algorithm:")
# start_time3 = time.time() * 1000000000
# print(s.KMP(pat, protein_seq))
# end_time3 = time.time() * 1000000000
#
# print("KMP takes ", start_time3 - end_time3, " to run")