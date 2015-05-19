import gzip
from xml.parsers import expat

class PDBMap(object):
    """Parse PDBML file for mappings between indexes of auth_seq_num
    + pdb_ind_code (residue index in PDB) and seq_id (starts with 1
    and includes unobserved residues).

    Note PDB indexes are in string, as it needs to combine with insertion
    code to unique and generally does not imply icreasing order.
    """
    def __init__(self, pdb):
        self.parser = expat.ParserCreate()
        self.parser.StartElementHandler = self._start_handler
        self.parser.EndElementHandler = self._end_handler
        self.parser.CharacterDataHandler = self._data_handler
        self.pdb = pdb.lower()
        self.pdbml_path = '/usr2/pdb/data/structures/divided/XML-noatom/' + \
                '%s/%s-noatom.xml.gz' % (self.pdb[1:3], self.pdb)
        self._tag = None
        self._chain = None
        self._seq_idx = None
        self._pdb_idx = None
        self._inscode = ''
        self._isordered = False
        self.mapping = {}
        with gzip.open(self.pdbml_path, 'rb') as fp:
            self.parser.ParseFile(fp)

    def _start_handler(self, name, attrs):
        if name == 'PDBx:pdbx_poly_seq_scheme':
            self._chain = attrs['asym_id']
            self._seq_idx = attrs['seq_id']
            if self._chain not in self.mapping:
                self.mapping[self._chain] = {'seq': [], 'pdb': []}
        elif self._chain and name == 'PDBx:auth_seq_num':
            self._isordered = True
            self._tag = 'auth_seq_num'
        elif self._chain and name == 'PDBx:pdb_ins_code':
            self._tag = 'pdb_ins_code'

    def _data_handler(self, data):
        if self._chain and self._tag == 'auth_seq_num':
            self._pdb_idx = data
        elif self._chain and self._tag == 'PDBx:pdb_ins_code':
            self._inscode = data

    def _end_handler(self, name):
        if name == 'PDBx:pdbx_poly_seq_scheme':
            if self._isordered:
                self.mapping[self._chain]['seq'].append(int(self._seq_idx))
                self.mapping[self._chain]['pdb'].append(self._pdb_idx +
                        self._inscode)
            self._isordered = False
            self._seq_id = None
            self._pdb_idx = None
            self._chain = None
            self._inscode = ''
        self._tag = None

    def get_mapping(self, chain):
        return self.mapping[chain]

# vim: ts=4 expandtab sw=4 sts=4 tw=78
