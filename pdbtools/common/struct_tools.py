from Bio.PDB import Model, Structure, PDBParser, PDBIO, parse_pdb_header, MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio import SeqIO
from Bio.Seq import Seq

import logging

from . import definitions


class Molecule(Model.Model):
    """Molecule is a wrapper around the Bio.PDB.Model object.
       It is effectivley the Model object but with added functionality.
       This functionality includes wrapped reader/writers,
       and Header parsing. Several other classes inherit from this."""
    def __init__(self, id, full_id, parent, child_list, child_dict, xtra,
                 serial_num = None, header = {}, seqres = None,
                 input_type = None):

        #Set Entity level properties
        self._id = id
        self.full_id = full_id
        self.parent = parent
        self.child_list = child_list
        self.child_dict = child_dict
        self.xtra = xtra

        #Set Model level properties
        self.level = "M"
        if serial_num is None:
            self.serial_num = id
        else:
            self.serial_num = serial_num

        #Additions
        self.header = header

        #Seqres should not have waters, and alts will ne converted or x
        self.seqres_seq = seqres

    @property
    def model_seq(self):
        """Fetches the sequence via the modeled atom records and returns as a dictionary"""
        seq_dict = {}
        for chain in self.get_chains():
            resns = [res.get_resname() for res in self[chain.id].get_residues()]
            resns = [definitions.three_to_one(res) for res in resns]
            chain_seq = Seq("".join(resns))
            seq_dict[chain.id] = chain_seq
        return seq_dict

    @property
    def model_resis(self):
        """Fetches the residue index via the modeled atom records and returns as a dictionary"""
        seq_dict = {}
        for chain in self.get_chains():
            resis = [res.id for res in self[chain.id].get_residues()]
            seq_dict[chain.id] = resis
        return seq_dict
            
    @classmethod
    def from_pdb(cls, file, id = "PDB_File"):
        '''Create a protein object from a '''
        #Parse the structure with the normal parser
        parser = PDBParser()
        rec_structure = parser.get_structure(id, file)
        if len(list(rec_structure.get_models())) > 1:
            logging.warning("More than one model (state) in pdb file -\
                            Using first one only")
        #Get the model of the structure
        rec_model = rec_structure[0]

        #Set the objects variables to those of the model (clone it)
        model_params = cls._get_model_params(rec_model)

        #Get the header information
        header = parse_pdb_header(file)
        seqres = cls.get_seqres(file)
        return cls(**model_params, header = header, seqres = seqres, input_type = "PDB")

    @classmethod
    def from_cif(cls, file, id = "CIF_File"):
        parser = MMCIFParser()
        structure = parser.get_structure(id, file)
        if len(list(structure.get_models())) > 1:
            logging.warning("More than one model (state) in pdb file -\
                            Using first one only")
        #Get the model of the structure
        rec_model = structure[0]

        #Set the objects variables to those of the model (clone it)
        model_params = cls._get_model_params(rec_model)

        #Get the header information
        header = MMCIF2Dict(file)
        seqres = cls.get_seqres(file, type = "cif")
        return cls(**model_params, header = header, seqres = seqres, input_type = "PDB")

    def renumber(self, start = 1, chain_restart = True, preserve_gaps = False):

        pass

    def renumber_to_seqres(self, start = 1):
        '''Renumpers the protein residues to match the seqres sequence'''
        pass

    #TODO
    @classmethod
    def from_structure(cls, struct):
        pass

    def write_pdb(self, filename):
        """Writes a PDB file from a biopython structure object"""
        io=PDBIO()
        io.set_structure(self)
        io.save(filename)

    def write_cif(self, filename):
        """Writes a PDB file from a biopython structure object"""
        io=MMCIFIO()
        io.set_structure(self)
        io.save(filename)

    @staticmethod
    def get_seqres(file, type = "pdb"):
        """Fetches the seqres record of a all chains in a structure file"""

        if type == "pdb":
            format = "pdb-seqres"
        elif type == "cif":
            format = "cif-seqres"
        else:
            raise ValueError("type: '{}' is not valid")
        seqres_dict = {}
        try:
            for record in SeqIO.parse(file, format):
                seqres_dict[record.annotations["chain"]] = Seq("".join(record.seq))
            return seqres_dict

        except Exception as e:
            return None

    @staticmethod
    def _get_model_params(model):
        """Takes a biopython model object and
           returns a dict with all the properties of that
           object. Used internally by several parsers"""
        id = model.id
        full_id = model.full_id
        parent = model.parent
        child_list = model.child_list
        child_dict = model.child_dict
        xtra = model.xtra
        serial_num = model.serial_num

        params = {"id":id,
                 "full_id":full_id,
                 "parent":parent,
                 "child_list":child_list,
                 "child_dict":child_dict,
                 "xtra":xtra,
                 "serial_num":serial_num
                 }
        return params

    def rename_nonstandard():
        pass

    def remove_nonstandard():
        pass

    def rebuild_nonstandard():
        pass

    def remove_hets():
        
        pass

    def renumber(self, numbering = "continuous", start = 1, chains = None):
        """Renumbers the residues, removing any insertion codes in the process
        args:
            numbering ['continuous', 'restart'] Determines whether the chains
                      should all restart numbering at 1, or continue from the
                      renumbering of the previous chain (default = continuous)
            start - Number to start renumbering on (default = int(1))
            chains - [None, str] Which chains to renumber. No argument will
                     renumber all chains. string arguents (eg. "A", "AC", will
                     only renumber those chains) (default = None)
        """
        if chains == None:
            chains = [chain.id for chain in self.get_chains()]
        else:
            chains = list(chains)

        # Avoid clobbering existing resids by adding "_temp"
        for chain in chains:
            for res in self[chain]:
                res.id = (res.id[0] + "_temp", res.id[1], res.id[2])

        #Then go through and renumber
        for chain in chains:
            for idx, res in enumerate(self[chain], start = start):
                res.id = (res.id[0].rstrip("_temp"),idx," ")

            if numbering == "continuous":
                start = start + len(self[chain])
            elif numbering == "restart":
                start = 1
            else:
                logging.warning("numbering value {} invalid".format(numbering))
                start = start + len(self[chain])



def rmsd(prot1, prot2):
    """Calcualtes the RMSD between two protein objects"""
    pass

def super(prot1, prot2, chain_map = None):
    """Superimposes two protein objects. First will be """
    pass

def cluster(*prot, cutoff = 1.0, method = "BB"):
    """Clusters protein objects to a certian cutoff in angstroms"""
    pass

def cluster_n(*prot, number = 3.,  method = "BB"):
    """Clusters protein objects into a certian number of clusters"""
    pass
