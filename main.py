from pdbtools.common import struct_tools


mymodel = struct_tools.Molecule.from_pdb("3WSO.pdb")
print(mymodel.seqres_seq["A"])