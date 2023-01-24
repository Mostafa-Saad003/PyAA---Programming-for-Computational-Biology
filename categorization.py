import characterization as c

def categorize_by_distance(protein_seq, list_of_sequences_1, list_of_sequences_2):

    # This main function serves to categorize a protein into one of two groups
    # ... by finding the difference between the composition of each amino acid in protein_seq and the average composition of the same amino acid in list_of_sequences_1 and list_of_sequences_2
    # Those differences are then summed, once for each group
    # If the distances difference is larger than the 6, normal composition is used
    # If the distances difference is less than 6, dipeptide composition is used to try and increase the divergence

    aa_composition = c.AA_Composition(protein_seq)
    average_composition_one = c.Avg_Composition(list_of_sequences_1)
    average_composition_two = c.Avg_Composition(list_of_sequences_2)

    category_1_distance = 0
    category_2_distance = 0

    for aa in aa_composition:
        category_1_distance += abs(aa_composition[aa] - average_composition_one[aa])
        category_2_distance += abs(aa_composition[aa] - average_composition_two[aa])

    if (abs((category_1_distance - category_2_distance)) > 6):
        if category_1_distance < category_2_distance:
            return "This protein belongs to category 1"
        else:
            return "This protein belongs to category 2"
    else:
        aa_pair_composition = c.pair_composition(protein_seq)
        average_pair_composition_one = c.average_pair_composition(list_of_sequences_1)
        average_pair_composition_two = c.average_pair_composition(list_of_sequences_2)

        category_1_distance = 0
        category_2_distance = 0

        for aa in aa_pair_composition:
            category_1_distance += abs(aa_pair_composition[aa] - average_pair_composition_one[aa])
            category_2_distance += abs(aa_pair_composition[aa] - average_pair_composition_two[aa])

        if category_1_distance < category_2_distance:
            return "This protein belongs to category 1"
        else:
            return "This protein belongs to category 2"

def folded_or_unfolded(protein_seq):
    # Two versions of code are provided below, both are inspired by the following paper: https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7
    # We coded the first version based on the explanation provided by the paper; however, its results are inaccurate.
    # ... This may be because we misinterpreted the calculation steps described in the paper or because the paper may have left out some details.
    # We coded the second version by trial and error by combining the steps mentioned in the paper into one equation.
    # ... This second method consistently produces accurate results, not just in terms of the final categorization but in terms of the mean hydrophobicity values themselves.

    # Kyte-Doolittle hydrophobicity scale
    aa_hydrophobicity = {'I': 4.50, 'V': 4.20, 'L': 3.80, 'F': 2.80, 'C': 2.50, 'M': 1.90, 'A': 1.80, 'G': -0.40, 'T': -0.70, 'S': -0.80,
                         'W': -0.90,'Y': -1.30, 'P': -1.60, 'H': -3.20, 'E': -3.50, 'N': -3.50, 'Q': -3.50, 'D': -3.50, 'K': -3.90, 'R': -4.50}

    # Version 1 of calculating the mean hydrophobicity values:
    # total_hydrophobicity = 0
    # To get total hydrophobicity, the window chosen by the paper was equal to 5 - this means that the average hydrophobicity for every 5 amino acids was calculated. These averages were summed to get the total hydrophobicity and divided by the protein length to get the required mean hydrophobicity.
    # half_window = 5//2
    # for i in range(half_window, len(protein_seq)-half_window, 5):
    #     window_sum = 0
    #     for j in range(-half_window, half_window+1):
    # The paper normalized the amino acid hydrophobicity values to 0 and 1 when getting their average over the window. They did not specify how they normalized; however, since the Kyte-Doolittle scale extends from -4.5 to 4.5, normalization was done by dividing hydrophobicity by 4.5 if hydrophobicity > 0 and -4.5 if hydrophobicity < 0.
    #         if aa_hydrophobicity[protein_seq[i+j]] < 0:
    #             normalized_aa = aa_hydrophobicity[protein_seq[i+j]]/-4.5
    #             window_sum += normalized_aa
    #         else:
    #             normalized_aa = aa_hydrophobicity[protein_seq[i+j]]/4.5
    #             window_sum += normalized_aa
    #     window_avg = window_sum/5
    #     total_hydrophobicity += window_sum
    # mean_hydrophobicity = total_hydrophobicity/len(protein_seq)

    # Version 2 (Brute-Force) of calculating the mean hydrophobicity values:
    mean_hydrophobicity = 2*(aa_hydrophobicity[protein_seq[0]] + aa_hydrophobicity[protein_seq[1]] + aa_hydrophobicity[protein_seq[-1]] + aa_hydrophobicity[protein_seq[-2]] + ((len(protein_seq)-4)/5))/(len(protein_seq))

    # Compares the mean hydrophobicity value calculated with the experimentally-derived hydrophobicity cutoff values of natively folded and unfolded proteins
    if (mean_hydrophobicity >= 0.34 and mean_hydrophobicity <= 0.44):
        return print("Protein is natively unfolded")
    elif (mean_hydrophobicity >= 0.45 and mean_hydrophobicity <= 0.51):
        return print("Protein is natively folded")
    else:
        return print('Inconclusive')
