from simphony.core.cuba import CUBA
from simphony.cuds.lattice import LatticeNode

from tables import *

import numpy
import copy


class FileLattice:
    """FileLattice object to use in file lattices.

    Parameters
    ----------
    file : CudsFile
        reference to PyTables file where the lattice is or will be located
    name : str
        name of the lattice
    type: string
        Bravais lattice type (should agree with the base_vect below).
    base_vect: D x float
        defines a Bravais lattice (an alternative for primitive vectors).
    size: D x size
        number of lattice nodes (in the direction of each axis).
    origin: D x float

    Returns
    ----------
    FileLattice
        The lattice newly added to the file or existing in the file.

    """

    class _ColDescr(IsDescription):
        pass

    def __init__(self, file, name, type, base_vect, size, origin):
        # Construct default FileLattice
        self._file = file._file
        self._group = file._file.root.lattice
        self._name = name
        self._type = type
        self._base_vect = base_vect
        self._size = size
        self._origin = origin
        self._maskname = name+'_bitmask'
        self._mask_size = 8

        # Construct default rowclass based on lattice metadata
        self._construct_default_rowclass()

        if self._name not in self._group:
            # Create default table to hold lattice data
            self._create_table()
            self._table.attrs.type = self._type
            self._table.attrs.base_vect = self._base_vect
            self._table.attrs.size = self._size
            self._table.attrs.origin = self._origin

        else:
            # Return reference to existing table
            self._table = self._file.root.lattice._f_getChild(self._name)
            self._type = self._table.attrs.type
            self._base_vect = self._table.attrs.base_vect
            self._size = self._table.attrs.size
            self._origin = self._table.attrs.origin
            # Update column description class.
            data = {}
            for key in self._table.coldescrs:
                if not key == 'NAME':
                    data[eval('CUBA.'+key)] = self._table[:][key][0]
            if not data == {}:
                self._update_coldescr(data, data.keys())

    def _construct_default_rowclass(self):
        """ Construct _ColDescr class with one column, so that
        PyTables table for the lattice can be generated. Empty table
        cannot be created.
        TODO: Find a more compact (and still meaningful) way of
        generating a dense table.
        """
        self._ColDescr = type("IsDescription", (IsDescription,),
                              {'NAME': Int32Col(pos=0,
                               shape=len(self._base_vect))})

    def _update_coldescr(self, data, keywords):
        """ Updates _ColDescr class to hold datatypes in DataContainers
        Parameters:
            data : dictionary
                Dictionary of CUBA-keywords with corresponding data
            keywords : list
                List of keywords that are added to _ColDescr class
        """
        # Keywords must be CUBA keywords. PyTables columns are named
        # after CUBA keywords where prefix 'CUBA.' is removed from the name.
        # TODO: Could use one predefined dictionary instead of the following
        # data type comparison
        # CD = self._ColDescr
        for key in keywords:
            # ND data has a different type than 1D
            if any([type(data[key]) == list,
                    type(data[key]) == numpy.ndarray]):
                if type(data[key][0]) == str:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: StringCol(16, dflt="",
                              pos=len(self._ColDescr.columns),
                              shape=len(data[key]))})
                elif type(data[key][0]) == long:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Int64Col(
                              pos=len(self._ColDescr.columns),
                              shape=len(data[key]))})
                elif type(data[key][0]) == float or numpy.float64:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Float64Col(pos=len(
                              self._ColDescr.columns),
                              shape=len(data[key]))})
                elif type(data[key][0]) == int or numpy.int32:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Int32Col(
                              pos=len(self._ColDescr.columns), shape=len(
                              data[key]))})
            else:
                if type(data[key]) == str:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: StringCol(16, dflt="",
                              pos=len(self._ColDescr.columns))})
                elif type(data[key]) == long:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Int64Col(pos=len(
                             self._ColDescr.columns))})
                elif type(data[key]) == int or numpy.int32:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Int32Col(pos=len(
                              self._ColDescr.columns))})
                elif type(data[key]) == float or numpy.float64:
                    self._ColDescr = type("_ColDescr", (self._ColDescr,),
                             {str(key)[5:]: Float64Col(pos=len(
                              self._ColDescr.columns))})

    def _index(self, n):
        """ Return ND index of row n, where N is the dimension of the lattice
        Parameters:
            n : int
                Row number
        Returns:
            I : tuple of int16
                ND tuple of lattice indices corresponding to row n
        """
        S = self._size
        Q = numpy.zeros(len(S),numpy.int64)
        I = numpy.zeros(len(S),numpy.int16)

        Q[len(S)-1] = 1
        for i in range(-len(S)+2,1):
            Q[-i] = Q[-i+1]*S[-i+1]

        I[0] = int(n/Q[0])
        for i in range(1,len(S)):
            I[i] = int((n-numpy.dot(I,Q))/Q[i])
        return tuple(I)

    def _get_row_index(self, id):
        """ Get row number from table for node with indices 'id'
        Parameters:
            id : tuple
                Indices of the LatticeNode
        Returns:
            n : int
                Row number
        """
        if len(self._size) == 3:
            n = numpy.dot(id, [self._size[1] * self._size[2],
                               self._size[2], 1])
            return n
        else:
            n = numpy.dot(id, [self._size[1], 1])
            return n

    def _create_table(self):
        """ Create default table and bitmask """
        F = self._file
        filters = Filters(complevel=5, complib='zlib')
        self._table = F.create_table(self._group, self._name,
                                     self._ColDescr, filters=filters)
        for n in xrange(numpy.prod(self._size)):
            row = self._table.row
            row['NAME'] = self._index(n)
            row.append()
        self._table.flush()
        # Max 8 data containers by default
        atom = UInt8Atom()
        self._bitmask = F.create_carray(self._group, self._maskname,
                                        atom, self._size, filters=filters)

    def _update_table(self, lattice_node, keywords):
        """ Checks if _ColDescr class already contains columns
        for DataContainers defined for lattice_node and adds them if
        necessary. Bitmask is updated as well.

        Parameters:
            lattice_node : LatticeNode
                reference to LatticeNode object
            keywords : list
                list of keywords that are updated to PyTables table

        """
        # Check that keywords have their dedicated columns in the table
        new_data_found = False
        for key in keywords:
            if str(key)[5:] not in self._table.coldescrs:
                new_data_found = True

        if new_data_found:
            # Construct a new table that can hold new data.
            # Column naming corresponds to CUBA keywords.
            # Because PyTables doesn't support '.' character
            # in column names 'CUBA.' is removed from the name string.
            new_keywords = set(copy.deepcopy(keywords))
            for key in self._table.coldescrs.keys():
                new_keywords = new_keywords.difference(
                    set([eval('CUBA.' + key)]))
            new_keywords = new_keywords.difference(set([CUBA(0)]))
            self._update_coldescr(lattice_node.data, new_keywords)
            filters = Filters(complevel=5, complib='zlib')
            newtable = self._file.create_table(self._group, "tmp",
                                               self._ColDescr,
                                               filters=filters)

            # Copy data from old table to new table
            for row in self._table.iterrows():
                filenode = newtable.row
                for key in self._table.coldescrs:
                    filenode[key] = copy.deepcopy(row[key])
                filenode.append()
            newtable.flush()

            # Copy attributes, must be done separately.
            newtable.attrs.type = copy.deepcopy(self._table.attrs.type)
            newtable.attrs.base_vect = \
                copy.deepcopy(self._table.attrs.base_vect)
            newtable.attrs.size = copy.deepcopy(self._table.attrs.size)
            newtable.attrs.origin = copy.deepcopy(self._table.attrs.origin)

            # Delete old table and rename new table
            self._table._f_remove()
            self._table = newtable
            self._file.rename_node(self._group, self._name, "tmp")

            # Check if larger bitmask has to be generated
            if len(self._ColDescr.columns.keys()) > self._mask_size:
                self._mask_size = 2*self._mask_size
                atom = eval('UInt'+str(self._mask_size)+'Atom()')
                filters = Filters(complevel=5, complib='zlib')
                self._newbitmask = \
                    self._file.create_carray(self._group,
                                             self._maskname+'_new', atom,
                                             self._size, filters=filters)
                self._newbitmask[:] = \
                    numpy.array(self._bitmask, dtype='uint'+str(
                                self._mask_size))

                self._bitmask._f_remove()
                self._bitmask = self._newbitmask
                self._file.root.lattice._v_leaves[
                    self._maskname+'_new'].rename(self._maskname)

    def _update_bitmask(self, id, name):
        """ Writes a bit 1 to bitmask to flag a defined DataContainer for
        lattice_node at position 'id'

        Parameters:
            id : tuple
                LatticeNode coordinates
            name : str
                Name of the PyTables column that is flagged
        """
        self._bitmask[id] = self._bitmask[id] |\
            pow(2, self._ColDescr.columns[name]._v_pos)

    def update_node(self, node, keywords=[]):
        """ Updates DataContainer contents for LatticeNode node

        Parameters:
            node : LatticeNode
                reference to LatticeNode object
            keywords : list, optional
                List of CUBA-keywords that are updated. If empty list
                then all DC's are from the node
        """
        # Find correct row for lattice_node
        n = self._get_row_index(node.id)

        # Check that table has all the needed columns
        if keywords == []:
            # Add everything
            self._update_table(node, node.data.keys())
            for row in self._table.iterrows(n, n+1):
                if len(node.data) is not 0:
                    for key in node.data:
                        name = str(key)[5:]
                        row[name] = copy.deepcopy(node.data[key])
                        self._update_bitmask(node.id, name)
                    row.update()
        else:
            # Add only selected keywords
            self._update_table(node, keywords)
            for row in self._table.iterrows(n, n+1):
                if len(node.data) is not 0:
                    for key in keywords:
                        name = str(key)[5:]
                        row[name] = copy.deepcopy(node.data[key])
                        self._update_bitmask(node.id, name)
                    row.update()

    def get_node(self, id):
        """Get a copy of the node corresponding to the given id.

        Parameters:
        -----------
        id: tuple of D x int (node index coordinate)

        Returns:
        -----------
        A reference to a LatticeNode object
        """
        n = self._get_row_index(id)
        tuple_id = tuple(id)
        dcs = {}
        for key in self._table.coldescrs.keys():
            dcs[eval('CUBA.'+key)] = self._table[n][key]
        del dcs[CUBA.NAME]
        return LatticeNode(tuple_id, dcs)

    def iter_nodes(self, ids=None):
        """Get an iterator over the LatticeNodes described by the ids.

        Parameters:
        -----------
        ids: iterable set of D x int (node index coordinates)

        Returns:
        -----------
        A generator for LatticeNode objects
        """
        if ids is None:
            # Iterate only over nodes that have some DataContainers
            # defined
            for id, val in numpy.ndenumerate(self._bitmask):
                if not val == 0:
                    yield self.get_node(id)
        else:
            for id in ids:
                yield self.get_node(id)

    def get_coordinate(self, id):
        """Get coordinate of the given index coordinate.

        Parameters:
        -----------
        id: D x int (node index coordinate)

        Returns:
        -----------
        D x float
        """
        return self._origin + self._base_vect*numpy.array(id)
