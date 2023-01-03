import requests
import numpy

def side_chain_charge(aa, pH=7):
    if aa == "A" or aa == "G" or aa == "V" or aa == "L" or aa == "I" or aa == "M" or aa == "P" or aa == "F" or aa == "W" or aa == "S" or aa == "T" or aa == "Q" or aa == "N":
        return 0
    elif aa == "C":
        if pH < 8.4:
            return 0
        elif pH == 8.4:
            return -0.5
        else:
            return -1
    elif aa == "Y":
        if pH < 10.5:
            return 0
        elif pH == 10.5:
            return -0.5
        else:
            return -1

    elif aa == "D" or aa == "E":
        if pH < 4:
            return 0
        elif pH == 4:
            return -0.5
        else:
            return -1

    elif aa == "K":
        if pH < 10.5:
            return 1
        elif pH == 10.5:
            return 0.5
        else:
            return 0

    elif aa == "R":
        if pH < 12.5:
            return 1
        elif pH == 12.5:
            return 0.5
        else:
            return 0

    elif aa == "H":
        if pH < 6:
            return 1
        elif pH == 6:
            return 0.5
        else:
            return 0
    else:
        return 0

def adjustTotal(num, protein_seq, pH):
    if protein_seq[-1] == "G" or protein_seq[-1] == "L" or protein_seq[-1] == "I" or protein_seq[-1] == "D" or protein_seq[-1] == "V":
        if pH > 9.6:
            num = num - 1
        elif pH == 9.6:
            num = num - 0.5
    elif protein_seq[-1] == "A" or protein_seq == "E":
        if pH > 9.7:
            num = num - 1
        elif pH == 9.7:
            num = num - 0.5
    elif protein_seq[-1] == "M" or protein_seq == "S" or protein_seq == "H":
        if pH > 9.2:
            num = num - 1
        elif pH == 9.2:
            num = num - 0.5
    elif protein_seq[-1] == "P":
        if pH > 10.6:
            num = num - 1
        elif pH == 10.6:
            num = num - 0.5
    elif protein_seq[-1] == "F" or protein_seq[-1] == "Q" or protein_seq[-1] == "T" or protein_seq[-1] == "Y":
        if pH > 9.1:
            num = num - 1
        elif pH == 9.1:
            num = num - 0.5
    elif protein_seq[-1] == "W":
        if pH > 9.4:
            num = num - 1
        elif pH == 9.4:
            num = num - 0.5
    elif protein_seq[-1] == "N":
        if pH > 8.8:
            num = num - 1
        elif pH == 8.8:
            num = num - 0.5
    elif protein_seq[-1] == "C":
        if pH > 8.2:
            num = num - 1
        elif pH == 8.2:
            num = num - 0.5
    elif protein_seq[-1] == "K" or protein_seq[-1] == "R":
        if pH > 9:
            num = num - 1
        elif pH == 9:
            num = num - 0.5
    return num

def parse_sequences(file):
    inFile = open(file, 'r')
    seqList = []
    currentSeq = ''
    for line in inFile:
        if line[0] == ">":
            if currentSeq != '':
                seqList.append(currentSeq)
            currentSeq = ''
        else:
            currentSeq += line.rstrip()
    seqList.append(currentSeq)

    return seqList

def fetch_proteins(URL, Name, num=1):
    with requests.get(URL, stream=True) as request:
        request.raise_for_status()
        n = Name + ".fasta"
        with open(n, 'wb') as f:
            for chunk in request.iter_content(chunk_size=2 ** 20):
                f.write(chunk)
    sequencesList = parse_sequences(n)
    if num == 1:
        sequencesList = "".join(sequencesList)
        return sequencesList
    else:
        return sequencesList

def pair_combinations():
    amino_acids = "ACDEFGHIKLMNPQRSTVWYU"
    combinations = {}
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            combination = amino_acids[i] + amino_acids[j]
            combinations[combination] = 0
    return combinations

pair_combs = pair_combinations()

def count_duplicate_pair(protein_seq, duplicate):
    sum = 0
    two_mers = len(protein_seq) - 2 + 1
    for i in range(two_mers):
        kmer = protein_seq[i: i+2]
        if kmer == duplicate:
            sum += 1
    return sum

def pair_occurence(protein_seq):
    pair_occurence = pair_combs
    for pair in pair_occurence:
        if pair[0] != pair[1]:
            pair_occurence[pair] = protein_seq.count(pair)
        else:
            pair_occurence[pair] = count_duplicate_pair(protein_seq, pair)
    return pair_occurence

def pair_composition(protein_seq):
    composition = {}
    for pair in pair_occurence(protein_seq):
        composition[pair] = round((pair_occurence(protein_seq)[pair] * 100) / (len(protein_seq) - 1), 2)
    return composition

def average_pair_composition(list_of_sequences):
    average_pair_composition = pair_combs
    # This is to increment the pairs
    for sequence in list_of_sequences:
        composition = pair_composition(sequence)
        for pair in composition:
            average_pair_composition[pair] += round(composition[pair], 2)

    # This is to find the average value of each pair
    for pair in average_pair_composition:
        average_pair_composition[pair] = round(average_pair_composition[pair] / len(list_of_sequences),2)

    return average_pair_composition

def isoelectric_point(protein_seq):
    min_positive = 100000
    min_negative = -100000
    pH_min_positive = 7
    pH_min_negative = 7
    zero_charge_list = []
    counter = 0
    for ph in numpy.arange(2.0, 12.6, 0.1):
        if total_charge(protein_seq, ph) == 0:
            counter += 1
            zero_charge_list.append(ph)
        elif total_charge(protein_seq, ph) > 0:
            if total_charge(protein_seq, ph) < min_positive:
                min_positive = total_charge(protein_seq, ph)
                pH_min_positive = ph
        else:
            if total_charge(protein_seq, ph) > min_negative:
                min_negative = total_charge(protein_seq, ph)
                pH_min_negative = ph

    if counter:
        return zero_charge_list
    else:
        return [pH_min_positive, pH_min_negative]

def average_isoelectric_point(list_of_sequences):
    summ = 0
    for seq in list_of_sequences:
        iso = sum(isoelectric_point(seq)) / len(isoelectric_point(seq))
        summ += iso

    return round(summ / len(list_of_sequences), 1)

def total_charge(protein_seq, pH=7):
    total = 0
    for aa in protein_seq:
        total += side_chain_charge(aa, pH)
    if pH < 2:
        total = total + 1
    elif pH == 2:
        total = total + 0.5
    else:
        total = adjustTotal(total, protein_seq, pH)
    return total

def average_total_charge(list_of_sequences, pH = 7):
    sum = 0
    for seq in list_of_sequences:
        sum += total_charge(seq, pH)

    return round(sum / len(list_of_sequences),1)

def num_negative_charge(protein_seq):
    return protein_seq.count("D") + protein_seq.count("E")

def average_num_negative_charge(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_negative_charge(seq)

    return round(sum / len(list_of_sequences),0)

def num_positive_charge(protein_seq):
    return protein_seq.count("R") + protein_seq.count("H") + protein_seq.count("K")

def average_num_positive_charge(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_positive_charge(seq)

    return round(sum / len(list_of_sequences),0)

def num_hydrophobic(protein_seq):
    return protein_seq.count("G") + protein_seq.count("A") + protein_seq.count("V") + protein_seq.count("L") + protein_seq.count("I") + protein_seq.count("P") + protein_seq.count("F") + protein_seq.count("M") + protein_seq.count("W")

def average_num_hydrophobic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_hydrophobic(seq)

    return round(sum / len(list_of_sequences),0)

def num_hydrophilic(protein_seq):
    return protein_seq.count("H") + protein_seq.count("E") + protein_seq.count("D") + protein_seq.count("N") + protein_seq.count("Q") + protein_seq.count("K") + protein_seq.count("R")

def average_num_hydrophilic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_hydrophilic(seq)

    return round(sum / len(list_of_sequences),0)

def num_aliphatic(protein_seq):
    return protein_seq.count("G") + protein_seq.count("A") + protein_seq.count("V") + protein_seq.count("L") + protein_seq.count("I") + protein_seq.count("P")

def average_num_aliphatic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_aliphatic(seq)

    return round(sum / len(list_of_sequences),0)

def num_aromatic(protein_seq):
    return protein_seq.count("F") + protein_seq.count("Y") + protein_seq.count("W")

def average_num_aromatic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_aromatic(seq)

    return round(sum / len(list_of_sequences),0)

def average_protein_length(list_of_sequences):
    sum = 0

    for seq in list_of_sequences:
        sum += len(seq)

    return round(sum / len(list_of_sequences),0)

def percent_negative_charge(protein_seq):
    return round(num_negative_charge(protein_seq) / len(protein_seq),2)

def average_percent_negative_charge(list_of_sequences):
    return round(average_num_negative_charge(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_positive_charge(protein_seq):
    return round(num_positive_charge(protein_seq) / len(protein_seq),2)

def average_percent_positive_charge(list_of_sequences):
    return round(average_num_positive_charge(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_hydrophobic(protein_seq):
    return round(num_hydrophobic(protein_seq) / len(protein_seq),2)

def average_percent_hydrophobic(list_of_sequences):
    return round(average_num_hydrophobic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_hydrophilic(protein_seq):
    return round(num_hydrophilic(protein_seq) / len(protein_seq),2)

def average_percent_hydrophilic(list_of_sequences):
    return round(average_num_hydrophilic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_aliphatic(protein_seq):
    return round(num_aliphatic(protein_seq) / len(protein_seq),2)

def average_percent_aliphatic(list_of_sequences):
    return round(average_num_aliphatic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_aromatic(protein_seq):
    return round(num_aromatic(protein_seq) / len(protein_seq),2)

def average_percent_aromatic(list_of_sequences):
    return round(average_num_aromatic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def export(input, name, file_name, num=1):
    exported_file = open(file_name, "w")
    amino_acids = "ACDEFGHIKLMNPQRSTVWYU"
    if num == 1:
        exported_file.writelines("Protein name: " + name + "\n")
        exported_file.writelines("Protein sequence: " + input + "\n")
        exported_file.writelines("Protein length: " + str(len(input)) + "\n")
        exported_file.writelines("Protein molecular weight: " + str(Protein_Molecular_Weight(input)) + "\n")
        exported_file.writelines("Protein net charge at pH = 7: " + str(total_charge(input, 7)) + "\n")
        exported_file.writelines("Protein isoelectric point range: " + str(round(isoelectric_point(input)[0],1)) + " - " + str(round(isoelectric_point(input)[-1],1)) + " | Mean: " + str(round(sum(isoelectric_point(input)) / len(isoelectric_point(input)), 2)) + "\n")
        exported_file.writelines("Protein net hydrophobicity: " + str(get_total_hydrophobicity(input)) + "\n")
        exported_file.writelines("Number of positively charged residues: " + str(num_positive_charge(input))+ " and their percentage: " + str(percent_positive_charge(input)) + "\n")
        exported_file.writelines("Number of negatively charged residues: "+ str(num_negative_charge(input)) + " and their percentage: " + str(percent_negative_charge(input)) + "\n")
        exported_file.writelines("Number of hydrophobic residues: "+ str(num_hydrophobic(input)) + " and their percentage: " + str(percent_hydrophobic(input)) + "\n")
        exported_file.writelines("Number of hydrophilic residues: "+ str(num_hydrophilic(input))+ " and their percentage: " + str(percent_hydrophilic(input)) + "\n")
        exported_file.writelines("Number of aliphatic residues: " + str(num_aliphatic(input)) + " and their percentage: " + str(percent_aliphatic(input)) + "\n")
        exported_file.writelines("Number of aromatic residues: " + str(num_aromatic(input)) + " and their percentage: " + str(percent_aromatic(input)) + "\n")

        exported_file.writelines("\n")
        exported_file.writelines("AA" + "\t" + "Occurence" + "\t" + "Composition" + "\n")
        occurence = AA_Occurence(input)
        composition = AA_Composition(input)
        for i in range(20):
            a = amino_acids[i]
            exported_file.writelines(a + "\t" + str(occurence[a]) + "\t" + str(composition[a]) + "\n")
            exported_file.writelines("\n")
    else:
        exported_file.writelines("Protein list name: " + name + "\n")
        exported_file.writelines("Dataset protein count: " + str(len(input)) + "\n")
        exported_file.writelines("Average protein length: " + str(average_protein_length(input)) + "\n")
        exported_file.writelines("Average protein molecular weight: " + str(Avg_molecular_weight(input)) + "\n")
        exported_file.writelines("Average protein net charge at pH = 7: " + str(average_total_charge(input, 7)) + "\n")
        exported_file.writelines("Average protein isoelectric point: " + str(average_isoelectric_point(input)) + "\n")
        exported_file.writelines("Average protein hydrophobicity: " + str(get_average_hydrophobicity(input)) + "\n")
        exported_file.writelines("Average number of positively charged residues: " + str(average_num_positive_charge(input))+ " and their percentage: " + str(average_percent_positive_charge(input)) + "\n")
        exported_file.writelines("Average number of negatively charged residues: "+ str(average_num_negative_charge(input)) + " and their percentage: " + str(average_percent_negative_charge(input)) + "\n")
        exported_file.writelines("Average number of hydrophobic residues: "+ str(average_num_hydrophobic(input)) + " and their percentage: " + str(average_percent_hydrophobic(input)) + "\n")
        exported_file.writelines("Average number of hydrophilic residues: "+ str(average_num_hydrophilic(input))+ " and their percentage: " + str(average_percent_hydrophilic(input)) + "\n")
        exported_file.writelines("Average number of aliphatic residues: " + str(average_num_aliphatic(input)) + " and their percentage: " + str(average_percent_aliphatic(input)) + "\n")
        exported_file.writelines("Average number of aromatic residues: " + str(average_num_aromatic(input)) + " and their percentage: " + str(average_percent_aromatic(input)) + "\n")
        exported_file.writelines("\n")
        exported_file.writelines("AA" + "\t" + "Occurence" + "\t" + "Composition" + "\n")
        avg_occurence = Avg_Occurence(input)
        avg_composition = Avg_Composition(input)
        for i in range(20):
            a = amino_acids[i]
            exported_file.writelines(a + "\t" + str(avg_occurence[a]) + "\t" + str(avg_composition[a]) + "\n")
            exported_file.writelines("\n")

def AA_Occurence(protein_seq):
    aa_occurence_dic = {}
    list_of_aa = ["A", "G", "P", "V", "L", "I",
                  "M", "C", "F", "Y", "W", "H",
                  "K", "R", "Q", "N", "E", "D",
                  "S", "T", "U"]
    for residue in list_of_aa:
        aa_occurence_dic[residue] = protein_seq.count(residue)
    return aa_occurence_dic

def Avg_Occurence(list_of_sequences):
    list_of_aa= ["A", "G", "P", "V", "L", "I",
                 "M", "C", "F", "Y", "W", "H",
                 "K", "R", "Q", "N", "E", "D",
                 "S", "T", "U"]

    average_occurence = {}
    for letter in list_of_aa:
        average_occurence[letter] = 0
    for seq in list_of_sequences:
        aa_occurence = AA_Occurence(seq)
        for aa in aa_occurence:
            average_occurence[aa] += aa_occurence[aa]

    for aa in average_occurence:
        average_occurence[aa] = average_occurence[aa] / len(list_of_sequences)

    return average_occurence

def AA_Composition(protein_seq):
    aa_composition_dic = AA_Occurence(protein_seq)
    for aa in aa_composition_dic:
        aa_composition_dic[aa] = round((aa_composition_dic[aa]/len(protein_seq))*100,2)
    return aa_composition_dic

def Avg_Composition(list_of_sequences):
    list_of_aa = ["A", "G", "P", "V", "L", "I",
                 "M", "C", "F", "Y", "W", "H",
                 "K", "R", "Q", "N", "E", "D",
                 "S", "T", "U"]

    average_composition = {}
    for letter in list_of_aa:
        average_composition[letter] = 0
    for seq in list_of_sequences:
        aa_composition = AA_Composition(seq)
        for aa in aa_composition:
            average_composition[aa] += aa_composition[aa]

    for aa in average_composition:
        average_composition[aa] = round(average_composition[aa] / len(list_of_sequences),2)

    return average_composition

def Protein_Molecular_Weight(protein_seq):
### Amino acids molecular weight source: https://worldwide.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/#:~:text=The%20average%20molecular%20weight%20of,(kDa)%20is%201%2C000%20daltons
    aa_mass_dic = {"A": 89, "G": 75, "P": 115, "V": 117, "L": 131, "I": 131,
                  "M": 149, "C": 121, "F": 165, "Y": 181, "W": 204, "H": 155,
                  "K": 146, "R": 174, "Q": 146, "N": 132, "E": 147, "D": 133,
                  "S": 105, "T": 119, "U": 167}
    protein_molecular_weight = 0
    for aa in protein_seq:
        protein_molecular_weight += aa_mass_dic[aa]

    protein_molecular_weight = round(protein_molecular_weight - (18 * (len(protein_seq) - 1)),2)
    return protein_molecular_weight

def Avg_molecular_weight(list_of_sequences):
    summ = 0
    for seq in list_of_sequences:
        summ += Protein_Molecular_Weight(seq)
    return round(summ / len(list_of_sequences),2)

def categorize_by_distance(protein_seq, list_of_sequences_1, list_of_sequences_2):
    aa_composition = AA_Composition(protein_seq)
    average_composition_one = Avg_Composition(list_of_sequences_1)
    average_composition_two = Avg_Composition(list_of_sequences_2)

    category_1_distance = 0
    category_2_distance = 0

    for aa in aa_composition:
        category_1_distance += abs(aa_composition[aa] - average_composition_one[aa])
        category_2_distance += abs(aa_composition[aa] - average_composition_two[aa])

    print(category_1_distance)
    print(category_2_distance)

    if (abs((category_1_distance - category_2_distance)) > 6):
        if category_1_distance < category_2_distance:
            print(abs(category_1_distance - category_2_distance))
            return "This protein belongs to category 1"
        else:
            return "This protein belongs to category 2"
    else:
        aa_pair_composition = pair_composition(protein_seq)
        average_pair_composition_one = average_pair_composition(list_of_sequences_1)
        average_pair_composition_two = average_pair_composition(list_of_sequences_2)

        category_1_distance = 0
        category_2_distance = 0

        for aa in aa_pair_composition:
            category_1_distance += abs(aa_pair_composition[aa] - average_pair_composition_one[aa])
            category_2_distance += abs(aa_pair_composition[aa] - average_pair_composition_two[aa])

        if category_1_distance < category_2_distance:
            return "This protein belongs to category 1"
        else:
            return "This protein belongs to category 2"

def get_total_hydrophobicity(protein_seq):
    aa_hydrophobicity = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37, 'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
                         'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62, 'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76}
    total_hydrophobicity = 0
    for aa in protein_seq:
        total_hydrophobicity += aa_hydrophobicity[aa]
    return round(total_hydrophobicity,2)

def get_average_hydrophobicity(list_of_proteins):
    sum_hydrophobicity = 0
    for protein in list_of_proteins:
        sum_hydrophobicity += get_total_hydrophobicity(protein)
    avg_hydrophobicity = sum_hydrophobicity/len(list_of_proteins)
    return round(avg_hydrophobicity,2)

def plot_hydrophobicity_profile(protein_seq, smoothing_window=1):
    import matplotlib.pyplot as plt
    aa_hydrophobicity = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37, 'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
                         'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62, 'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76}
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
    #plt.ylim(-2, 2)
    plt.show()

def composition_plot(protein, num=1):
    import matplotlib.pyplot as plt
    plt.title("Protein Composition")
    if num == 1:
        composition = AA_Composition(protein)
    else:
        composition = Avg_Composition(protein)
    plt.xlabel("Amino Acid")
    plt.ylabel("Frequency (%)")
    amino_acid = list(composition.keys())
    frequency = list(composition.values())
    plt.bar(amino_acid, frequency, color ='orange', width = 0.4)
    for i in range(len(amino_acid)):
        plt.text(i, frequency[i], frequency[i], ha = 'center')

    plt.show()

def categorize_by_correlation(protein_seq, list_of_sequences1, list_of_sequences2):
    prot_composition = AA_Composition(protein_seq)
    category1 = AA_Composition(list_of_sequences1)
    category2 = AA_Composition(list_of_sequences2)
    numerator = 0
    denominator = 0
    sum1_sq = 0
    sum2_sq = 0
    for aa in prot_composition:
        sum1_sq += ((prot_composition[aa] - 5)**2)
        sum2_sq += ((category1[aa] - 5) ** 2)
        numerator += ((prot_composition[aa] - 5) * (category1[aa] - 5))
    denominator = ((sum1_sq * sum2_sq) ** 0.5)
    r1 = numerator / denominator

    numerator = 0
    denominator = 0
    sum1_sq = 0
    sum2_sq = 0
    for aa in prot_composition:
        sum1_sq += ((prot_composition[aa] - 5) ** 2)
        sum2_sq += ((category2[aa] - 5) ** 2)
        numerator += ((prot_composition[aa] - 5) * (category2[aa] - 5))
    denominator = ((sum1_sq * sum2_sq) ** 0.5)
    r2 = numerator / denominator

    if r1 > r2:
        return 'This protein belongs to category 1'
    elif r2 > r1:
        return 'This protein belongs to category 2'

#transmembrane_proteins = fetch_proteins("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%29%20AND%20%28proteins_with%3A56%29%20AND%20%28proteins_with%3A1%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29%20AND%20%28existence%3A1%29%20AND%20%28annotation_score%3A5%29", "Transmembrane", 1354)
#VDAC = fetch_proteins("https://rest.uniprot.org/uniprotkb/P21796.fasta", "VDAC", 1)
# dna_binding = fetch_proteins("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29%20AND%20%28model_organism%3A9606%29%20AND%20%28proteins_with%3A1%29%20AND%20%28proteins_with%3A23%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28annotation_score%3A5%29%20AND%20%28length%3A%5B201%20TO%20400%5D%29", "DNA-Binding Proteins", 58)
# #motif_proteins = fetch_proteins("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28proteins_with%3A11%29%20AND%20%28proteins_with%3A35%29", "Motif", 11)
# #tata_binding_protein = fetch_proteins("https://rest.uniprot.org/uniprotkb/P20226.fasta", "Tata-Binding Protein", 1)
#TFIID = fetch_proteins("https://rest.uniprot.org/uniprotkb/Q8IZX4.fasta", "TFIID", 1)
#receptor = fetch_proteins("https://rest.uniprot.org/uniprotkb/Q16620.fasta", "Receptor", 1)
#plot_hydrophobicity_profile(VDAC, 20)

tfiidd = "MRPGCDLLLRAAATVTAAIMSDSDSEEDSSGGGPFTLAGILFGNISGAGQLEGESVLDDECKKHLAGLGALGLGSLITELTANEELTGTGGALVNDEGWIRSTEDAVDYSDINEVAEDESQRHQQTMGSLQPLYHSDYDEDDYDADCEDIDCKLMPPPPPPPGPMKKDKDQDAITCVSESGEDIILPSIIAPSFLASEKVDFSSYSDSESEMGPQEATQAESEDGKLTLPLAGIMQHDATKLLPSVTELFPEFRPGKVLRFLHLFGPGKNVPSVWRSARRKRKKHRELIQEEQIQEVECSVESEVSQKSLWNYDYAPPPPPEQCLADDEITMMVPVESKFSQSTGDVDKVTDTKPRVAEWRYGPARLWYDMLGVSEDGSGFDYGFKLRKTQHEPVIKSRMMEEFRKLEESNGTDLLADENFLMVTQLHWEDSIIWDGEDIKHKGTKPQGASLAGWLPSIKTRNVMAYNVQQGFAPTLDDDKPWYSIFPIDNEDLVYGRWEDNIIWDAQAMPRLLEPPVLALDPNDENLILEIPDEKEEATSNSPSKESKKESSLKKSRILLGKTGVIREEPQQNMSQPEVKDPWNLSNDEYYFPKQQGLRGTFGGNIIQHSIPAMELWQPFFPTHMGPIKIRQFHRPPLKKYSFGALSQPGPHSVQPLLKHIKKKAKMREQERQASGGGELFFMRTPQDLTGKDGDLILAEYSEENGPLMMQVGMATKIKNYYKRKPGKDPGAPDCKYGETVYCHTSPFLGSLHPGQLLQALENNLFRAPVYLHKMPETDFLIIRTRQGYYIRELVDIFVVGQQCPLFEVPGPNSRRANMHIRDFLQVFIYRLFWKSKDRPRRIRMEDIKKAFPSHSESSIRKRLKLCADFKRTGMDSNWWVLKSDFRLPTEEEIRAKVSPEQCCAYYSMIAAKQRLKDAGYGEKSFFAPEEENEEDFQMKIDDEVHAAPWNTTRAFIAAMKGKCLLEVTGVADPTGCGEGFSYVKIPNKPTQQKDDKEPQAVKKTVTGTDADLRRLSLKNAKQLLRKFGVPEEEIKKLSRWEVIDVVRTMSTEQAHSGEGPMSKFARGSRFSVAEHQERYKEECQRIFDLQNKVLSSTEVLSTDTDSISAEDSDFEEMGKNIENMLQNKKTSSQLSREWEEQERKELRRMLLVAGSAASGNNHRDDVTASMTSLKSSATGHCLKIYRTFRDEEGKEYVRCETVRKPAVIDAYVRIRTTKDEKFIQKFALFDEKHREEMRKERRRIQEQLRRLKRNQEKEKLKGPPEKKPKKMKERPDLKLKCGACGAIGHMRTNKFCPLYYQTNVPPSKPVAMTEEQEEELEKTVIHNDNEELIKVEGTKIVFGKQLIENVHEVRRKSLVLKFPKQQLPPKKKRRVGTTVHCDYLNIPHKSIHRRRTDPMVTLSSILESIINDMRDLPNTHPFHTPVNAKVVKDYYKIITRPMDLQTLRENVRKCLYPSREEFREHLELIVKNSATYNGPKHSLTQISQSMLDLCDEKLKEKEDKLARLEKAINPLLDDDDQVAFSFILDNIVTQKMMAVPDSWPFHHPVNKKFVPDYYKMIVNPVDLETIRKNISKHKYQSRESFLDDVNLILANSVKYNGPESQYTKTAQEIVNICYQTITEYDEHLTQLEKDICTAKEAALEEAELESLDPMTPGPYTSQPPDMYDTNTSLSTSRDASVFQDESNLSVLDISTATPEKQMCQGQGRLGEEDSDVDVEGYDDEEEDGKPKPPAPEGGDGDLADEEEGTVQQPEASVLYEDLLISEGEDDEEDAGSDEEGDNPFSAIQLSESGSDSDVGYGGIRPKQPFMLQHASGEHKDGHGK"
receptorr = "MSSWIRWHGPAMARLWGFCWLVVGFWRAAFACPTSCKCSASRIWCSDPSPGIVAFPRLEPNSVDPENITEIFIANQKRLEIINEDDVEAYVGLRNLTIVDSGLKFVAHKAFLKNSNLQHINFTRNKLTSLSRKHFRHLDLSELILVGNPFTCSCDIMWIKTLQEAKSSPDTQDLYCLNESSKNIPLANLQIPNCGLPSANLAAPNLTVEEGKSITLSCSVAGDPVPNMYWDVGNLVSKHMNETSHTQGSLRITNISSDDSGKQISCVAENLVGEDQDSVNLTVHFAPTITFLESPTSDHHWCIPFTVKGNPKPALQWFYNGAILNESKYICTKIHVTNHTEYHGCLQLDNPTHMNNGDYTLIAKNEYGKDEKQISAHFMGWPGIDDGANPNYPDVIYEDYGTAANDIGDTTNRSNEIPSTDVTDKTGREHLSVYAVVVIASVVGFCLLVMLFLLKLARHSKFGMKGPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIENPQYFGITNSQLKPDTFVQHIKRHNIVLKRELGEGAFGKVFLAECYNLCPEQDKILVAVKTLKDASDNARKDFHREAELLTNLQHEHIVKFYGVCVEGDPLIMVFEYMKHGDLNKFLRAHGPDAVLMAEGNPPTELTQSQMLHIAQQIAAGMVYLASQHFVHRDLATRNCLVGENLLVKIGDFGMSRDVYSTDYYRVGGHTMLPIRWMPPESIMYRKFTTESDVWSLGVVLWEIFTYGKQPWYQLSNNEVIECITQGRVLQRPRTCPQEVYELMLGCWQREPHMRKNIKGIHTLLQNLAKASPVYLDILG"
vdacc = "MAVPPTYADLGKSARDVFTKGYGFGLIKLDLKTKSENGLEFTSSGSANTETTKVTGSLETKYRWTEYGLTFTEKWNTDNTLGTEITVEDQLARGLKLTFDSSFSPNTGKKNAKIKTGYKREHINLGCDMDFDIAGPSIRGALVLGYEGWLAGYQMNFETAKSRVTQSNFAVGYKTDEFQLHTNVNDGTEFGGSIYQKVNKKLETAVNLAWTAGNSNTRFGIAAKYQIDPDACFSAKVNNSSLIGLGYTQTLKPGIKLTLSALLDGKNVNAGGHKLGLGLEFQA"

list_of_proteins = ["MAVPPTYADLGKSA", "MAVPPTYADLGKSAA", "MAVPPTYADLGKSAAA"]
export(list_of_proteins, "List", "list.txt", 3)
export(vdacc, "VDAC", "v.txt", 1)