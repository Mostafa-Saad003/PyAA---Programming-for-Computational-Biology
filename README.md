# BMS 321 - Programming For Computational Biology - PyAA 1.0

# What is PyAA?

This is a python package that helps molecular biologists, structural biologists, and protein bioinformaticians do comprehensive analyses on protein sequences. 

# How can the package help me?

The package consists of four primary modules: 1) characterization, 2) categorization, 3) structure_prediction, and 4) searching

1) characterization --> Analyzes the amino acid sequence of a protein or a list of proteins (average) --> i.e. total charge, isoelectric point, molecular weight, composition, etc.
2) categorization --> Categorizes a protein, either to one of any two protein groups, or based on the native folding status of the protein (folded vs. unfolded)
That is, VDAC for example, can be categorized twice --> Transmembrane vs. DNA-binding protein OR Folded vs. Unfolded
3) structure_prediction --> Predicts the secondary structure of the protein based on 3 methods: ramachandran plot of the protein's dihedral angles, Chou and Fasman's statistical propensity method, and the hydrophobicity plot of the protein
4) searching --> Searches with a pattern (motif of interest) in a protein sequence based on 3 algorithms: Naive, Boyer-Moore, and KMP

# How to use the package?

To use the package, simply import the moudles that contains the analysis of interest
i.e. if you want to export a txt file for protein characterization, "import characterization" and proceed with your desired analysis

To see a real example, check below "How to run a demo?"

# How to run a demo?

'Main.py' has a demo for every module of the four. To test each, just uncomment the related code.

Created By: Ahmed Hesham Saadawy, Lydia Reda Sidarous, Mostafa Saad Ibrahim
