#!/usr/bin/env python3
import pandas as pd
# import argparse
# import sys
# import os
# import itertools
# import math
from Bio.PDB import *
# from numpy import *
# from scipy.spatial import cKDTree
# from Bio import BiopythonWarning
# import warnings
# from Bio.Data.SCOPData import protein_letters_3to1
# from Bio.SeqUtils import seq1
# from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo
# from Bio import SeqIO
import logging

class capristat():
    protein_backbone = ["N", "CA", "C", "O"]
    nucleic_backbone = ["C1'", "C2'", "C3'", "C4'", "C5'", "O3'", "O4'", "O5'", "OP1", "OP2", "P"]
    total_backbone = protein_backbone + nucleic_backbone

    protein_resns = list("GALMFWKQESPVICYHRNDTX")
    nucleic_resns = list("gactu")

    def __init__(self, models, reference,rchains=None, mrchains=None, lchains=None, mlchains=None):
        #Input paramaters and data
        self.reference = reference
        self.models = models
        self.rchains = rchains
        self.mrchains = mrchains
        self.lchains = lchains
        self.mlchains = mlchains
        #Dataframe to store final metrics
        self.metrics = pd.DataFrame()
        #Dataframe to store contact info
        self.contacts = pd.DataFrame()

        #Refrence chain information and sequences
        self.refchains =  list(chain.get_id() for chain in self.reference.get_chains()) #List of chain ID's
        self.refseqs = {}
        for chain in self.refchains:
            self.refseqs[chain] = []
            for residue in self.refchains.get_residues():
                self.refseqs[chain].append(residue.)

    def map_chains(self, gap_open = -10, gap_extend = -0.5):
        '''Takes unaligned '''
        pass

    def get_clashes(self, cutoff = 3):
        pass
    def get_fnat(self):
        pass
    def get_irmsd(self):
        pass
    def get_lrmsd(self):
        pass
    def get_quality(self):
        pass

    def write_pymol(self):
        pass


    @classmethod
    def from_filepaths(cls, reference_path, *modelpaths):
        pass

    @classmethod
    def from_commandline(cls):
        pass

# if __name__ == "__main__":
#     capri.from_commandline()

parser = PDBParser(PERMISSIVE=True)
ref = parser.get_structure("ref", "1A2K_c_b.pdb")       #These are STRUCTURES
struct1 = parser.get_structure("model1", "1A2K_c_u.pdb")

refmodel = ref[0]       #These are the forst MODEL (frame/state)
model1 = struct1[0]

capristats = capristat(model1, reference = refmodel)
