from simphony.core.cuba import CUBA

from tables import *

import numpy
import copy


class FileLattice:

    class _LatticeDescription(IsDescription):
        pass

    def __init__(self, file, name='lattice', lattice=None):
        """Returns FileLattice object.

        Parameters
        ----------
        file : PyTables file
            reference to PyTables file where the lattice is or will be located
        name : str, optional
            name of the lattice, default name is 'lattice'
        lattice : Lattice, optional
            lattice to be added. If none is given,
            then an empty lattice called 'lattice' is added. If lattice
            with name 'name' already exists, then return a reference to
            that lattice.

        Returns
        ----------
        FileLattice
            The lattice newly added to the file or existing in the file.

        """

        # Construct default FileLattice
        self._file = file
        self._group = file.root.lattice
        self._name = name
        self._type = 'Cubic'
        self._base_vect = (1, 1, 1)
        self._size = (10, 10, 10)
        self._origin = (0, 0, 0)

        if lattice is not None:
            # Construct FileLattice based on Lattice
            self._name = lattice.name
            self._type = lattice.type
            self._base_vect = lattice.base_vect
            self._size = lattice.size
            self._origin = lattice.origin

        # Construct default rowclass based on lattice metadata
        self._construct_default_rowclass()

        if self._name not in self._group:
            # Create default table to hold lattice data
            self._create_lattice_table()
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
            # Update rowclass: ugly implementation, but will have to do atm
            data = {}
            for key in self._table.coldescrs:
                if not key == 'NAME':
                    data[eval('CUBA.'+key)] = self._table[:][key][0]
            if not data == {}:
                self._update_rowclass(data, data.keys())

    def _construct_default_rowclass(self):
        self._LatticeDescription = type("IsDescription", (IsDescription,),
                                        {'NAME': Int32Col(pos=0,
                                         shape=len(self._base_vect))})

    def _update_rowclass(self, data, keywords):
        # Keywords must be CUBA keywords. PyTables columns are named
        # after CUBA keywords where prefix 'CUBA.' is removed from the name.
        for key in keywords:
            # ND data has a different type than 1D
            if any([type(data[key]) == list,
                    type(data[key]) == numpy.ndarray]):
                LD = self._LatticeDescription
                if type(data[key][0]) == str:
                    LD = type("_LatticeDescription", (LD,),
                              {str(key)[5:]: StringCol(16, dflt="",
                               pos=len(LD.__dict__)+1,
                               shape=len(data[key]))})
                elif type(data[key][0]) == long:
                    LD = type("_LatticeDescription", (LD,),
                              {str(key)[5:]: Int64Col(
                               pos=len(LD.__dict__)+1,
                               shape=len(data[key]))})
                elif type(data[key][0]) == int or numpy.int32:
                    LD = type("_LatticeDescription", (LD,),
                              {str(key)[5:]: Int32Col(
                               pos=len(LD.__dict__)+1,
                               shape=len(data[key]))})
                elif type(data[key][0]) == float or numpy.float64:
                    LD = type("_LatticeDescription", (LD,),
                              {str(key)[5:]: Float64Col(
                               pos=len(LD.__dict__)+1,
                               shape=len(data[key]))})
            else:
                if type(data[key]) == str:
                    LD = type("_LatticeDescription", (LD,), {str(key)[5:]:
                              StringCol(16, dflt="",
                              pos=len(LD.__dict__)+1)})
                elif type(data[key]) == long:
                    LD = type("_LatticeDescription", (LD,), {str(key)[5:]:
                              Int64Col(pos=len(LD.__dict__)+1)})
                elif type(data[key]) == int or numpy.int32:
                    LD = type("_LatticeDescription", (LD,), {str(key)[5:]:
                              Int32Col(pos=len(LD.__dict__)+1)})
                elif type(data[key]) == float or numpy.float64:
                    LD = type("_LatticeDescription", (LD,), {str(key)[5:]:
                              Float64Col(pos=len(LD.__dict__)+1)})

    def _create_lattice_table(self):
        self._table = self._file.create_table(self._group, self._name,
                                              self._LatticeDescription)

    def _update_table(self, lattice_node, keywords):
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
                new_keywords.difference(set([eval('CUBA.' + key)]))
            new_keywords.difference(set([CUBA(0)]))
            self._update_rowclass(lattice_node.data, new_keywords)
            newtable = self._file.create_table(self._group, "tmp",
                                               self._LatticeDescription)
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

    def update_node(self, lattice_node, keywords=[]):
        # First check that table has all the needed columns
        if keywords == []:
            # Add everything
            self._update_table(lattice_node, lattice_node.data.keys())
        else:
            self._update_table(lattice_node, keywords)

        # Find correct row for lattice_node
        updated = False
        for row in self._table:
            if tuple(row["NAME"]) == lattice_node.id:
                for key in lattice_node.data:
                    row[str(key)[5:]] = copy.deepcopy(lattice_node.data[key])
                # Fix this
                row.update()
                self._table.flush()
                break
        if not updated:
            self.add_node(lattice_node, keywords)

    def add_node(self, lattice_node, keywords=[]):
        # lattice_node.data must have keywords defined!
        # For performance reasons it is not
        # checked every time a node is added.

        # Create new node in the table. (use update node,
        # if node exists already)
        if not keywords == []:
            self._update_table(lattice_node, keywords)

        filenode = self._table.row
        filenode["NAME"] = copy.deepcopy(lattice_node.id)
        if not keywords == []:
            self._update_table(lattice_node, keywords)
            for key in keywords:
                filenode[str(key)[5:]] = \
                copy.deepcopy(lattice_node.data[key])
        filenode.append()
        self._table.flush()
