import logging
from . import definitions


class Sequence():
    '''Biopython Seq object with added Functionality'''
    def __init__():
        pass


def read_fasta(fasta_file):
    """Reads a single entry fasta file and returns a biopython
     sequence object"""
    seqs = list(SeqIO.parse(fasta_file, "fasta"))
    if len(seqs) > 1:
        logging.warning("More than one sequence in fasta file -\
                        Using first one only")
    seq = seqs[0] #.seq
    return seq

def read_multifasta(fasta_file):
    """Reads a multi entry fasta file and returns a biopython
     sequence object"""
    seqs = list(SeqIO.parse(fasta_file, "fasta"))
    return seqs

def get_sequence_string(seq_obj):
    """Converts a biopython sequence record object to a string"""
    sequence_string = "".join([resn for resn in seq_obj])
    return sequence_string

def write_fasta_string(sequence, id = "sequence", description = "fasta file",
                output = "sequence.fasta"):
    """Writes a fasta file from a sequence string"""
    print(type(sequence))
    fasta_record = SeqRecord(Seq(sequence), id=id, description=description)
    with open(output, "w") as output_handle:
        SeqIO.write(fasta_record, output_handle, "fasta")

def write_fasta(sequence,  output = "sequence.fasta"):
    """Writes a fasta file from a biopython sequence object"""
    with open(output, "w") as output_handle:
        SeqIO.write(sequence, output_handle, "fasta")


warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import parse_pdb_header, PDBParser, Polypeptide
from Bio.PDB.PDBIO import PDBIO
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import warnings
from Bio import BiopythonWarning


from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")


def get_pdbseq_polypeptide(model, target_chain):
    """Extract chain sequence via biopython polypeptide builder"""
    ppb=Polypeptide.PPBuilder()
    pp = ppb.build_peptides(model[target_chain])
    seq = Seq("".join([str(frag.get_sequence()) for frag in pp]))
    return seq

def renumber_chain(model, target_chain, aln_target, aln_pdb):
    """Renumbers a pdb chain to match the numbering an a target alignment. Will
    Delete PDB residues not present in the alignment.
    model           : Biopython model
    target_chain    : Which chain to renumber
    aln_target      : The alignment to match to - likely the SEQRES / FASTA file
    aln_pdb         : The PDB chain sequence to align to the target. Will be renumbered
    """
    #Create a dictionary to match the {PDB seq index: Target seq index}
    pdb_to_seq = {}
    pdb_idx = 0
    seq_idx = 0
    for idx, (seqaln, pdbaln) in enumerate(zip(aln_target, aln_pdb)):
        if seqaln != "-": #IF seqres seq is present in this position
            seq_idx += 1
            deposit_seq = seq_idx
        else:
            deposit_seq = None
        if pdbaln != "-": #IF pdb sequence is present in this position
            pdb_idx += 1
            deposit_pdb = pdb_idx
            pdb_to_seq[pdb_idx] = deposit_seq
        else:
            deposit_pdb = None

    #Avoid clobbering existing resids
    for res in model[target_chain]:
        res.id = (res.id[0] + "_temp", res.id[1], res.id[2])

    #Renumber the reamining ones to match their aligned fasta index
    remove_ids = []
    for idx, res in enumerate(model[target_chain], start = 1):
        new_id = pdb_to_seq[idx]
        if new_id == None: #Not in the Seqres records - Remove
            remove_ids.append(res.id)
        else:
            res.id = (res.id[0].rstrip("_temp"),new_id," ") #Set to the true ID

    for res_id in remove_ids:
        model[target_chain].detach_child(res_id)

def global_align(seq_1, seq_2):
    """Gets global alignment between two sequences, returns aligned sequences
       and score """
    alignments = pairwise2.align.globalds(seq_1, seq_2, blosum62, -10, -0.5,
                                          penalize_end_gaps=False)
    if len(alignments) == 0:
        print("Alignment Failed!")
        return None
    best_align = sorted(alignments, key = lambda x: x.score, reverse = True)[0]
    return best_align.seqA, best_align.seqB, best_align.score

def read_fasta(fasta_file):
    """Reads a single entry fasta file and returns a biopython
     sequence object"""
    seqs = list(SeqIO.parse(fasta_file, "fasta"))
    if len(seqs) > 1:
        logging.warning("More than one sequence in fasta file -\
                        Using first one only")
    seq = seqs[0] #.seq
    return seq

def get_sequence_string(seq_obj):
    """Converts a biopython sequence record object to a string"""
    sequence_string = "".join([resn for resn in seq_obj])
    return sequence_string

def write_fasta_string(sequence, id = "sequence", description = "fasta file",
                       output = "sequence.fasta"):
    """Writes a fasta file from a sequence string"""
    print(type(sequence))
    fasta_record = SeqRecord(Seq(sequence), id=id, description=description)
    with open(output, "w") as output_handle:
        SeqIO.write(fasta_record, output_handle, "fasta")

def write_fasta(sequence,  output = "sequence.fasta"):
    """Writes a fasta file from a biopython sequence object"""
    with open(output, "w") as output_handle:
        SeqIO.write(sequence, output_handle, "fasta")
