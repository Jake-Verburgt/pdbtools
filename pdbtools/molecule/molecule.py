from Bio.PDB import Model, Structure, PDBParser, PDBIO, parse_pdb_header
from Bio.PDB.StructureBuilder import StructureBuilder



prot_three_to_one = {
    'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
    'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
    'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
    'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

prot_one_to_three = {item[1]: item[0] for item in three_to_one.items()}

class Molecule(Model.Model):
    """Protein is a wrapper around the Bio.PDB.Model object.
       It is effectivley a biopython Model object but with added functionality.
       This functionality includes wrapped reader/writers,
       and Header parsing. Several other classes inherit from this."""


   protein_backbone = ["N", "CA", "C", "O"]
   nucleic_backbone = ["C1'", "C2'", "C3'", "C4'", "C5'", "O3'", "O4'", "O5'",
                       "OP1", "OP2", "P"]
   total_backbone = protein_backbone + nucleic_backbone

    def __init__(self, id, full_id, parent, child_list, child_dict, xtra,
                 serial_num = None, header = {}
                 ):

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


    @classmethod
    def from_structure(cls, struct):
        pass

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


    def write_pdb(self, filename):
        """Writes a PDB file from a biopython structure object"""
        io=PDBIO()
        io.set_structure(self)
        io.save(filename)

    @classmethod
    def from_pdb(cls, file, id = "PDB_File"):
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

        return cls(**model_params, header = header)

        def get_nucleic_chains():
            pass


        def get_protein_chains():
            pass

        def

class Complex(Molecule):
    """Extension of Protein class for handling complex-specific tasks
    'Complex' does not mean multi-chain, but rather a ligand/receptor model"""


    def from_molecule()
        pass
