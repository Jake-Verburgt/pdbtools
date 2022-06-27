
def three_to_one(code, standardize = True):
    '''Coverts a 3 letter code to 1 letter code, with support for waters and hets'''
    if code in water_three_to_one.keys():
        return water_three_to_one[code]
    if code in protein_three_to_one.keys():
        return protein_three_to_one[code]
    elif code in nucleic_three_to_one.keys():
        return nucleic_three_to_one[code]
    elif code in protein_substitutions.keys() and standardize:
        return protein_three_to_one[protein_substitutions[code]]
    elif code in nucleic_substitutions.keys() and standardize:
        return nucleic_three_to_one[nucleic_substitutions[code]]
    else:
        return "X"


# Protein definitions
protein_three_to_one = {
    'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
    'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
    'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
    'GLY': 'G', 'PRO': 'P', 'CYS': 'C'
}

protein_one_to_three = {item[1]: item[0] for item in protein_three_to_one.items()}
protein_res_3 = list(protein_three_to_one.keys())
protein_res_1 = list(protein_three_to_one.values())


protein_substitutions = {
    '2AS': 'ASP', '3AH': 'HIS', '5HP': 'GLU', 'ACL': 'ARG', 'AGM': 'ARG', 'AIB': 'ALA', 'ALM': 'ALA', 'ALO': 'THR', 'ALY': 'LYS', 'ARM': 'ARG',
    'ASA': 'ASP', 'ASB': 'ASP', 'ASK': 'ASP', 'ASL': 'ASP', 'ASQ': 'ASP', 'AYA': 'ALA', 'BCS': 'CYS', 'BHD': 'ASP', 'BMT': 'THR', 'BNN': 'ALA',
    'BUC': 'CYS', 'BUG': 'LEU', 'C5C': 'CYS', 'C6C': 'CYS', 'CAS': 'CYS', 'CCS': 'CYS', 'CEA': 'CYS', 'CGU': 'GLU', 'CHG': 'ALA', 'CLE': 'LEU', 'CME': 'CYS',
    'CSD': 'ALA', 'CSO': 'CYS', 'CSP': 'CYS', 'CSS': 'CYS', 'CSW': 'CYS', 'CSX': 'CYS', 'CXM': 'MET', 'CY1': 'CYS', 'CY3': 'CYS', 'CYG': 'CYS',
    'CYM': 'CYS', 'CYQ': 'CYS', 'DAH': 'PHE', 'DAL': 'ALA', 'DAR': 'ARG', 'DAS': 'ASP', 'DCY': 'CYS', 'DGL': 'GLU', 'DGN': 'GLN', 'DHA': 'ALA',
    'DHI': 'HIS', 'DIL': 'ILE', 'DIV': 'VAL', 'DLE': 'LEU', 'DLY': 'LYS', 'DNP': 'ALA', 'DPN': 'PHE', 'DPR': 'PRO', 'DSN': 'SER', 'DSP': 'ASP',
    'DTH': 'THR', 'DTR': 'TRP', 'DTY': 'TYR', 'DVA': 'VAL', 'EFC': 'CYS', 'FLA': 'ALA', 'FME': 'MET', 'GGL': 'GLU', 'GL3': 'GLY', 'GLZ': 'GLY',
    'GMA': 'GLU', 'GSC': 'GLY', 'HAC': 'ALA', 'HAR': 'ARG', 'HIC': 'HIS', 'HIP': 'HIS', 'HMR': 'ARG', 'HPQ': 'PHE', 'HTR': 'TRP', 'HYP': 'PRO',
    'IAS': 'ASP', 'IIL': 'ILE', 'IYR': 'TYR', 'KCX': 'LYS', 'LLP': 'LYS', 'LLY': 'LYS', 'LTR': 'TRP', 'LYM': 'LYS', 'LYZ': 'LYS', 'MAA': 'ALA', 'MEN': 'ASN',
    'MHS': 'HIS', 'MIS': 'SER', 'MLE': 'LEU', 'MPQ': 'GLY', 'MSA': 'GLY', 'MSE': 'MET', 'MVA': 'VAL', 'NEM': 'HIS', 'NEP': 'HIS', 'NLE': 'LEU',
    'NLN': 'LEU', 'NLP': 'LEU', 'NMC': 'GLY', 'OAS': 'SER', 'OCS': 'CYS', 'OMT': 'MET', 'PAQ': 'TYR', 'PCA': 'GLU', 'PEC': 'CYS', 'PHI': 'PHE',
    'PHL': 'PHE', 'PR3': 'CYS', 'PRR': 'ALA', 'PTR': 'TYR', 'PYX': 'CYS', 'SAC': 'SER', 'SAR': 'GLY', 'SCH': 'CYS', 'SCS': 'CYS', 'SCY': 'CYS',
    'SEL': 'SER', 'SEP': 'SER', 'SET': 'SER', 'SHC': 'CYS', 'SHR': 'LYS', 'SMC': 'CYS', 'SOC': 'CYS', 'STY': 'TYR', 'SVA': 'SER', 'TIH': 'ALA',
    'TPL': 'TRP', 'TPO': 'THR', 'TPQ': 'ALA', 'TRG': 'LYS', 'TRO': 'TRP', 'TYB': 'TYR', 'TYI': 'TYR', 'TYQ': 'TYR', 'TYS': 'TYR', 'TYY': 'TYR',
    'HSD': 'HIS', 'HSE': 'HIS', 'HIE': 'HIS', 'CYX': 'CYS', 'MP8': 'PRO', 'HIP': 'HIS', 'CYZ': 'CYS'
}


protein_backbone_atoms = ["N", "CA", "C", "O"]
protein_sidechain_atoms= {
    "ALA":["CB"],
    "ARG":["CB", "CG","CD", "NE", "CZ", "NH1", "NH2"],
    "ASN":["CB", "CG", "OD1", "ND2"],
    "ASP":["CB", "CG", "OD1", "OD2"],
    "CYS":["CB", "SG"],
    "GLU":["CB", "CG", "CD", "OE1", "OE2"],
    "GLN":["CB", "CG", "CD", "OE1", "NE2"],
    "GLY":[],
    "HIS":["CB", "CG", "ND1", "CE1", "NE2", "CD2"],
    "ILE":["CB", "CG1", "CG2", "CD1"],
    "LEU":["CB", "CG", "CD1", "CD2"],
    "LYS":["CB", "CG", "CD", "CE", "NZ"],
    "MET":["CB", "CG", "SD", "CE"],
    "PHE":["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO":["CB", "CG", "CD"],
    "SER":["CB", "OG"],
    "THR":["CB", "OG1", "CG2"],
    "TRP":["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR":["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL":["CB", "CG1", "CG2"]
}


protein_backbone_bonds = [["N","CA"], ["CA", "C"], ["C","O"]]
protein_sidechain_bonds = {
    "ALA":[["CA","CB"]],
    "ARG":[["CA","CB"], ["CB","CG"], ["CG","CD"], ["CD", "NE"], ["NE","CZ"], ["CZ","NH1"], ["CZ","NH2"]],
    "ASN":[["CA","CB"], ["CB","CG"], ["CG","OD1"], ["CG","ND2"]],
    "ASP":[["CA","CB"], ["CB","CG"], ["CG","OD1"], ["CG","OD2"]],
    "CYS":[["CA","CB"], ["CB","SG"]],
    "GLU":[["CA","CB"], ["CB","CG"], ["CG","CD"], ["CD","OE1"], ["CD","OE2"]],
    "GLN":[["CA","CB"], ["CB","CG"], ["CG","CD"], ["CD","OE1"], ["CD","NE2"]],
    "GLY":[],
    "HIS":[["CA","CB"], ["CB","CG"], ["CG","ND1"], ["ND1","CE1"], ["CE1","NE2"], ["NE2","CD2"], ["CD2", "CG"]],
    "ILE":[["CA","CB"], ["CB","CG1"], ["CG1","CD1"], ["CB","CG2"]],
    "LEU":[["CA","CB"], ["CB", "CG"], ["CG", "CD1"], ["CG","CD2"]],
    "LYS":[["CA","CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "CE"], ["CE","NZ"]],
    "MET":[["CA","CB"], ["CB", "CG"], ["CG", "SD"], ["SD","CE"]],
    "PHE":[["CA","CB"], ["CB", "CG"], ["CG", "CD1"], ["CG1","CE1"],["CE1", "CZ"], ["CG", "CD2"], ["CG2","CE2"], ["CE2", "CZ"]],
    "PRO":[["CA","CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "N"]],
    "SER":[["CA","CB"], ["CB", "OG"]],
    "THR":[["CA","CB"], ["CB", "OG1"], ["CB", "OG2"]],
    "TRP":[["CA","CB"], ["CB", "CG"], ["CG", "CD1"], ["CD1", "NE1"],["NE1", "CE2"], ["CE2", "CZ2"], ["CZ2", "CH2"], ["CH2", "CZ3"], ["CG", "CD2"], ["CD2", "CE2"],["CD2", "CE3"], ["CE3", "CZ3"]],
    "TYR":[["CA","CB"], ["CB", "CG"], ["CG", "CD1"], ["CG1","CE1"],["CE1", "CZ"], ["CG", "CD2"], ["CG2","CE2"], ["CE2", "CZ"], ["CZ", "OH"]],
    "VAL":[["CA","CB"], ["CB", "CG1"], ["CB", "CG2"]]
}


protein_terminals = ["OXT"]
protein_terminal_alts = {"OT1": "O", 
                         "OT2": "OXT"}


#Mapping CHARMM protein atom names back to standard
charmm_protein = {}


# Nucleic definitions
nucleic_three_to_one = {"DA":"A", "A":"A",
                        "DG":"G", "G":"G",
                        "DC":"C", "C":"C",
                        "DT":"T", "T":"T",
                        "DU":"U", "U":"U"
}

nucleic_substitutions = {}

nucleic_backbone_atoms = ["P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*"
                          "O3*", "C2*", "O2*", "C1*"]

nucleic_sidechain_atoms = {"DA":["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2",
                                 "N3", "C4"],
                           "DG":["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2",
                                 "N2", "N3", "C4"],
                           "DC":["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],

                           "DT":["N1", "C2", "O2", "N3", "C4", "O4", "C5"," C5M",
                                 "C6"],
                           "U" :["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
                           }

#Mapping CHARMM nucleic atom names back to standard
charmm_nucleic = {}


# Water Definitions
water_three_to_one = {"HOH":"O", "OH2":"O", "H2O":"O","WAT":"O"}

