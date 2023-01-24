import math
import numpy as np
import matplotlib.pyplot as plt


# Dictionaries that store Chou and Fasman's Parameters
alpha_propensity_dictionary = {"A": [1, 1.42], "C": [-1, 0.70], "D": [0.5, 1.01], "E": [1, 1.39], "F": [1, 1.13],
                                   "G": [-1, 0.57], "H": [0.5, 1.00], "I": [0.5, 1.08], "K": [1, 1.14], "L": [1, 1.41],
                                   "M": [1, 1.45], "N": [-1, 0.67], "P": [-1, 0.57], "Q": [1, 1.11], "R": [0, 0.98],
                                   "S": [0, 0.77], "T": [0, 0.83], "W": [0.5, 1.08], "V": [0.5, 1.06], "Y": [-1, 0.69]}
beta_propensity_dictionary = {"A": [0, 0.83], "C": [1, 1.19], "D": [-1, 0.54], "E": [1, 1.17], "F": [1, 1.38],
                                  "G": [0, 0.75], "H": [0, 0.87], "I": [1, 1.60], "K": [0, 0.74], "L": [1, 1.30],
                                  "M": [0.5, 1.05], "N": [0, 0.89], "P": [-1, 0.55], "Q": [0.5, 1.10], "R": [0, 0.93],
                                  "S": [0, 0.75], "T": [1, 1.19], "W": [1, 1.37], "V": [1, 1.70], "Y": [1, 1.47]}

def get_secondary_structures_by_chou_and_fasman(seq):

    # Main function to get alpha helices and beta sheets using Chou and Fasman's propensity method

    dict_of_alpha = {}
    dict_of_beta = {}
    helix = ""
    beta = ""
    j_alpha = 0
    j_beta = 0

    # I am interested in finding all the 6-mers and 5-mers, but since 6-mer is larger I use it as the limit of the loop

    kmer_len = len(seq) - 6 + 1
    for i in range(kmer_len):
        score_alpha = 0 # to keep track of relative amino acid alpha helix score
        absolute_score_alpha = 0 # to keep track of absolute amino acid alpha helix score
        score_beta = 0 # to keep track of relative amino acid beta sheet score
        absolute_score_beta = 0 # to keep track of absolute amino acid beta sheet score
        sixmer = seq[i: i+6] # the 6-mer
        fivemer = seq[i: i+5] # the 5-mer
        j_alpha = i # to store the start of the 6-mer, representing alpha helix start
        j_beta = i # to store the start of the 5-mer, representing beta sheet start

        # Calculating the total relative scores of the 6-mer (score_alpha_ and 5-mer (score_beta)
        for amino_acid in sixmer:
            score_alpha += alpha_propensity_dictionary[amino_acid][0]
        for amino_acid in fivemer:
            score_beta += beta_propensity_dictionary[amino_acid][0]

        # This condition means the segment is likely dominantly an alpha helix, and not a beta sheet
        if score_alpha >= 4 and score_beta < 3:
            counter = 0 # to trace any additions on the original 6-mers
            helix = sixmer # originally, the helix is the 6-mer
            for j in range(j_alpha+6, kmer_len): # we loop over other 4-mers to keep increasing the helix segment, until it's no longer a helix
                # No longer a helix means that the absolute_score of a 4-mer is < 4
                temp_score = 0 # to keep track of absolute alpha score
                temp_score_beta = 0 # to keep track of absolute beta score
                fourmer = seq[j: j+4] # the 4-mer
                for amino_acid in fourmer:
                    temp_score += alpha_propensity_dictionary[amino_acid][1]
                    temp_score_beta += beta_propensity_dictionary[amino_acid][1]

                # This condition checks that the segment is strongly alpha helix (temp_score >= 4), but that it is also stronger than a beta sheet (temp_score > temp_score_beta)
                # The checking of counter == 0 is to check whether this is the first time we're elongating the segment or not
                # If it is the first time, we add the entire 4-mer
                if (temp_score >= 4) and (temp_score > temp_score_beta) and counter == 0:
                    counter += 1
                    helix = sixmer + fourmer

                # However, if it's not the first time of elongation, we only add the last residue, because those are 4-mers, meaning the second 4-mer already has 3 residues from the 4-mer preceding it
                elif (temp_score >= 4) and (temp_score > temp_score_beta) and counter > 0:
                    helix += fourmer[-1]

                # If the segment to be elongated is not a strong helix former, we append the last saved helix segment to a dictionary of alpha helices, where the key is the segment and the value is the start/end of the segment
                else:
                    start = j_alpha + 1 # start of the segment --> + 1 because we begin the loop from the very beginning from 0
                    end = j_alpha + len(helix) # to know where the segment ends, just add the start to the length of the segment :'D
                    dict_of_alpha[helix] = (start, end) # this is the key/value pair of the dictionary
                    break # the loop is broken so that the outer loop is accessed

        # The upcoming conditions repeat the same logic but according to whether the region is a helix former or a beta former

        elif score_beta >= 3 and score_alpha < 4:
            counter = 0
            beta = fivemer
            for j in range(j_beta+5, kmer_len):
                temp_score = 0
                temp_score_alpha = 0
                trimer = seq[j: j + 3]
                for amino_acid in trimer:
                    temp_score += beta_propensity_dictionary[amino_acid][1]
                    temp_score_alpha += alpha_propensity_dictionary[amino_acid][1]
                if (temp_score >= 3) and (temp_score > temp_score_alpha) and counter == 0:
                    counter += 1
                    beta = fivemer + trimer
                elif (temp_score >= 3) and (temp_score > temp_score_alpha) and counter > 0:
                    beta += trimer[-1]
                else:
                    start = j_beta + 1
                    end = j_beta + len(beta)
                    dict_of_beta[beta] = (start, end)
                    break

        # if the region is both a strong alpha helix former and a beta sheet former, we check the average of the 6-mer and 5-mer scores to see which is larger
        # If the 6-mer average "absolute" score is larger, it means it's alpha-helix
        # If the 5-mer average "absolute" score is larger, it means it's beta sheet

        elif (score_alpha >= 4 and score_beta >= 3) and ((absolute_score_alpha / 6) > (absolute_score_beta / 5)):
            counter = 0
            helix = sixmer
            for j in range(j_alpha+6, kmer_len):
                temp_score = 0
                temp_score_beta = 0
                fourmer = seq[j: j + 4]
                for amino_acid in fourmer:
                    temp_score += alpha_propensity_dictionary[amino_acid][1]
                    temp_score_beta += beta_propensity_dictionary[amino_acid][1]
                if (temp_score >= 4) and (temp_score > temp_score_beta) and counter == 0:
                    counter += 1
                    helix = sixmer + fourmer
                elif (temp_score >= 4) and (temp_score > temp_score_beta) and counter > 0:
                    helix += fourmer[-1]
                else:
                    start = j_alpha + 1
                    end = j_alpha + len(helix)
                    dict_of_alpha[helix] = (start, end)
                    break
        elif (score_alpha >= 4 and score_beta >= 3) and ((absolute_score_alpha / 6) < (absolute_score_beta / 5)):
            counter = 0
            beta = fivemer
            for j in range(j_beta+5, kmer_len):
                temp_score = 0
                temp_score_alpha = 0
                trimer = seq[j: j + 3]
                for amino_acid in trimer:
                    temp_score += beta_propensity_dictionary[amino_acid][1]
                    temp_score_alpha += alpha_propensity_dictionary[amino_acid][1]
                if (temp_score >= 3) and (temp_score > temp_score_alpha) and counter == 0:
                    counter += 1
                    beta = fivemer + trimer
                elif (temp_score >= 3) and (temp_score > temp_score_alpha) and counter > 0:
                    beta += trimer[-1]
                else:
                    start = j_beta + 1
                    end = j_beta + len(beta)
                    dict_of_beta[beta] = (start, end)
                    break

        # This should have conditions for other secondary structures, such as coils and turns.
        # But for now, it will just proceed with the next loop without appending anything to the dictionary
        else:
            continue

    # Eventually, the output is a dictionary where the key is either the word "Alpha Helices" or "Beta Sheets"
    # And the value of each key is another dictionary
    # That other dictionary contains key/value pairs of "Segment: (start, end)"
    return {"Alpha Helices": dict_of_alpha, "Beta Sheets": dict_of_beta}

def modify_chou_and_fasman(seq):

    # Because the function above will yield many overlapping segments, i.e. 6-15 and 7-12 both are alpha-helices
    # We need to remove all the segments that are substrings of other longer segments
    # This function serves to do that: to remove all keys that are substrings of other keys

    d = get_secondary_structures_by_chou_and_fasman(seq)
    # Create a list of keys to delete
    keys_to_delete = []
    # Iterate over the keys in the dictionary
    for key1 in d["Alpha Helices"]:
        # Iterate over the keys again
        for key2 in d["Alpha Helices"]:
            # If key1 is a substring of key2 and key1 is not equal to key2
            if key1 in key2 and key1 != key2:
                # Add key1 to the list of keys to delete
                keys_to_delete.append(key1)
    for key1 in d["Beta Sheets"]:
        # Iterate over the keys again
        for key2 in d["Beta Sheets"]:
            # If key1 is a substring of key2 and key1 is not equal to key2
            if key1 in key2 and key1 != key2:
                # Add key1 to the list of keys to delete
                keys_to_delete.append(key1)
    # Delete the keys from the dictionary
    for key in keys_to_delete:
        if key in d["Alpha Helices"]:
            del d["Alpha Helices"][key]
        elif key in d["Beta Sheets"]:
            del d["Beta Sheets"][key]
        else:
            continue

    return d

def generate_DSSP_sequence(seq):

    # This function mainly iterates over the sequence and compares the alpha/beta propensity values of each amino acid
    # If the alpha-value > beta-value AND also it's a strong helix former ( > 1), it writes "H" over the same amino acid position
    # If the beta-value > alpha-value AND also it's a strong beta former ( > 1), it writes "E" over the same amino acid position
    # If none of those, we don't do anything, so eventually a "C" will be written --> Coil
    # Eventually, we have a string that looks like this --> "EEHHHHCCHEHEHHCCCCHEEEEEECCCHHHHHHHHE", etc...
    # This method of representing a protein sequence is called a DSSP sequence

    # Initialize the prediction array
    prediction = np.zeros(len(seq))
    # Iterate over the sequence and predict the secondary structure of each amino acid
    for i, aa in enumerate(seq):
        h = alpha_propensity_dictionary[aa][1]
        e = beta_propensity_dictionary[aa][1]
        if h > e and h >= 1:
            prediction[i] = 1
        elif h < e and e >= 1:
            prediction[i] = -1
    # Convert the prediction array to a string of characters
    prediction = ''.join(['H' if x == 1 else 'E' if x == -1 else 'C' for x in prediction])
    return prediction

def getCoordinates(pdb_file):

    # This function gives the xyz coordinates (as a list of tuples) of all the N, CA, and C atoms in a pdb file

    # Initialize a list to store the coordinates
    coordinates = []
    # Open the PDB file
    with open(pdb_file, 'r') as f:
        # Iterate over each line of the file
        for line in f:
            # Check if the line contains atom coordinate data
            if line.startswith('ATOM'):
                elements = line.split() # to split the line into columns
                if elements[2] == "N" or elements[2] == "CA" or elements[2] == "C": # Bec. I am only interested in those atoms
                    # Extract the x, y, and z coordinates from the line
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # Add the coordinates to the list
                    coordinates.append((x, y, z))
    return coordinates

def calculate_cos(v1, v2, v3):

    # This function calculates cos theta of the vectors v1 x v2 and v2 x v3
    # To know what v1, v2, and v3 will be, proceed with the usage of this function below

    # Take 3 vectors as inputs

    # Get new vectors, which come from the cross product of the original vectors
    v4 = np.cross(v1, v2)
    v5 = np.cross(v2, v3)

    # Calculate the dot product of the vectors
    dot_product = v4[0] * v5[0] + v4[1] * v5[1] + v4[2] * v5[2]

    # Calculate the magnitudes of the vectors
    v4_magnitude = math.sqrt(v4[0] ** 2 + v4[1] ** 2 + v4[2] ** 2)
    v5_magnitude = math.sqrt(v5[0] ** 2 + v5[1] ** 2 + v5[2] ** 2)

    # Calculate the angle between the vectors using the dot product and magnitudes
    cos = dot_product / (v4_magnitude * v5_magnitude)

    return cos

def calculate_sin(v1, v2, v3):

    # This function calculates sin theta of the vectors v1 x v2 and v2 x v3
    # This is because the dihedral angle is defined in the range (-pi, pi), so we need sin theta as well
    # We need sin theta so that we can use it with cos theta to get atan2(sin, cos) to know the angle with its approporiate sign
    # We could have used cos only, but we'll need further checks to confirm if it's positive or negative
    # To know what v1, v2, and v3 will be, proceed with the usage of this function below

    # Take 3 vectors as inputs

    # Get new vectors, which come from the cross product of the original vectors
    v4 = np.cross(v1, v2)
    v5 = np.cross(v2, v3)

    # Calculate the cross product of the new crossed vectors and dot it with v2
    cross_product = np.cross(v4, v5)
    dot_product = cross_product[0] * v2[0] + cross_product[1] * v2[1] + cross_product[2] * v2[2]

    # Calculate the magnitudes of the vectors
    v4_magnitude = math.sqrt(v4[0] ** 2 + v4[1] ** 2 + v4[2] ** 2)
    v5_magnitude = math.sqrt(v5[0] ** 2 + v5[1] ** 2 + v5[2] ** 2)
    v2_magnitude = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)

    # Calculate the angle between the vectors using the dot product and magnitudes
    sin = dot_product / (v4_magnitude * v5_magnitude * v2_magnitude)

    return sin

def get_residue_coordinates(coordinates, residue_number):

    # This function takes coodrinates as inputs (list of tuples) as well as a residue number
    # It outputs the coordinates of this particulr residue obly

    return coordinates[3*(residue_number-1): residue_number*3]

def subtract_vectors(v1,v2):

    # This function is to subtract two vectors from each other

    return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])

def calculate_phi_psi(sequence, coordinates):

    # The main function of calculating phi and psi

    angles = [] # Saving the angles here, as tuples that contain 2 elements: phi and psi
    for i in range(1, len(sequence) - 1): # Looping from 1 not 0 because the first residue does not have a phi --> does not have a previous residue
        # and till len(seq) - 1 because the last residue does not have a psi --> does not have a next
        # Calculate the phi angle for the current residue
        prev_res_coords = get_residue_coordinates(coordinates, i) # previous residue coordinates
        curr_res_coords = get_residue_coordinates(coordinates, i+1) # current residue coordinates
        next_res_coords = get_residue_coordinates(coordinates, i+2) # next residue coordinates

        N = curr_res_coords[0] # Current residue nitrogen
        C_prev = prev_res_coords[2] # Previous residue C carbonyl
        C_current = curr_res_coords[2] # Current residue C carbonyl
        Ca = curr_res_coords[1] # Current residue alpha Carbon
        N_next = next_res_coords[0] # Next residue nitrogen

        # For calculating phi
        a = subtract_vectors(N, C_prev)
        b = subtract_vectors(Ca, N)
        c = subtract_vectors(C_current, Ca)

        cos_phi = calculate_cos(a, b, c) # the vectors to which we calculate cos corresponding to phi
        sin_phi = calculate_sin(a, b, c) # the vectors to which we calculate sin corresponding to psi
        phi = math.atan2(sin_phi, cos_phi) * 57.2958 # because it's in radians

        # For calculating psi
        if i < len(sequence) - 2: # to check for the next out-of-index issue
            # By the same token as above, what differs is just the three vectors to which we calculate cos and sin for.
            # They differ as in different connectivities, because the psi bond connects differnt atoms than what the phi does.
            a_psi = subtract_vectors(Ca, N)
            b_psi = subtract_vectors(C_current, Ca)
            c_psi = subtract_vectors(N_next, C_current)
            cos_psi = calculate_cos(a_psi, b_psi, c_psi)
            sin_psi = calculate_sin(a_psi, b_psi, c_psi)
            psi = math.atan2(sin_psi, cos_psi) * 57.2958  # because it's in radians
        else:
            psi = None

        angles.append((phi, psi))

    return angles

def plot_ramachandran(list_of_phi_psi):

    # This function serves to graph the Ramachandran plot for the dihedral angles of a protein sequence
    # X axis has the range that the Phi angle values can take (-180, 180)
    # Y axis has the range that the Psi angle values can take (-180, 180)
    # The dihedral angles calculated for every amino acid are plotted as points
    # X-coordinate --> Phi angle, Y-coordinate --> Psi angle

    plt.title('Ramachandran Plot')
    plt.xlabel('Phi Φ Angle')
    plt.xlim(-180, 180)
    plt.ylabel('Psi Ψ Angle')
    plt.ylim(-180, 180)
    list_of_phi_psi.pop()
    plt.hlines(y=0, xmin= -180, xmax= 180)
    plt.vlines(x=0, ymin= -180, ymax= 180)
    for dihedral in list_of_phi_psi:
        # beta sheets were found to have dihedral angles corresponding to this range of values
        if (dihedral[0] <= -80 and dihedral[0] >= -150) and (dihedral[1] >= 70 and dihedral[1] <= 170):
            plt.scatter(dihedral[0], dihedral[1], color='red')
        # alpha helices were found to have dihedral angles corresponding to this range of values
        elif (dihedral[0] <= -50 and dihedral[0] >= -150) and (dihedral[1] >= -60 and dihedral[1] <= -10):
            plt.scatter(dihedral[0], dihedral[1], color='orange')
        else:
            plt.scatter(dihedral[0], dihedral[1], color='blue')
    plt.show()

def plot_hydrophobicity_profile(protein_seq, smoothing_window=1):

    # This function serves to plot a hydrophobicity profile for a protein sequence
    # X axis has the amino acids, sequentially arranged
    # Y axis has the hydrophobic value for that particular amino acid
    # These dots are then connected to form the hydrophobic profile

    import matplotlib.pyplot as plt
    aa_hydrophobicity = {'I': 4.50, 'V': 4.20, 'L': 3.80, 'F': 2.80, 'C': 2.50, 'M': 1.90, 'A': 1.80, 'G': -0.40, 'T': -0.70, 'S': -0.80,
                         'W': -0.90,'Y': -1.30, 'P': -1.60, 'H': -3.20, 'E': -3.50, 'N': -3.50, 'Q': -3.50, 'D': -3.50, 'K': -3.90, 'R': -4.50}
    plt.title('Hydrophobicity Profile')
    plt.xlabel('Amino acid sequence')

    x_values = []
    y_values = []
    if smoothing_window == 1:
        plt.ylabel('Hydrophobicity index')
        for aa in range(len(protein_seq)):
            x_values.append(aa)
            y_values.append(aa_hydrophobicity[protein_seq[aa]])
        plt.plot(x_values, y_values)

    # Smoothing window plots the average hydrophobicity value of the amino acids within the odd window

    elif smoothing_window > 1:
        for aa in range(len(protein_seq)):
            y_values.append(aa_hydrophobicity[protein_seq[aa]])
        half_window = smoothing_window//2
        plt.ylabel('Hydrophobicity index (window = %d)' % smoothing_window)
        y_average = []
        for i in range(half_window, len(protein_seq)-half_window):
            average = 0.0
            x_values.append(i)
            for j in range(-half_window, half_window+1):
                average += y_values[i+j]
            y_average.append(average/smoothing_window)
        plt.plot(x_values, y_average)
    plt.show()