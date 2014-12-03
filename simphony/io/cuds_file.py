""" File access to CUDS-hdf5 formatted files

This module provides access (read/write) to file which are
formated in the CUDS-hdf5 file format
"""

import copy

import tables

from simphony.io.file_particle_container import FileParticleContainer
from file_lattice import FileLattice


class CudsFile(object):
    """ Access to CUDS-hdf5 formatted files

    """
    def __init__(self, file):
        """

        Parameters
        ----------
        file : table.file
            file to be used

        """

        if not isinstance(file, tables.File):
            raise ValueError(
                "File should be a Pytable file")

        if file.mode is 'r':
            raise ValueError(
                "File should not be opened in read-only mode")

        self._file = file
        self._particle_containers = {}
        self._lattices = {}

    def valid(self):
        """Checks if file is valid (i.e. open)

        """
        return self._file is not None and self._file.isopen

    @classmethod
    def open(cls, filename, mode="a", title=''):
        """Returns a SimPhony file and returns an opened CudsFile

        Parameters
        ----------
        filename : str
            Name of file to be opened.


        mode: str
            The mode to open the file:
                * *'w'*: Write; a new file is created (an existing file
                  with the same name would be deleted).
                * *'a'*: Append; an existing file is opened for reading and
                  writing, and if the file does not exist it is created.


        title : str
            Title attribute of root node (only applies to a file which
              is being created

        """
        if mode not in ('a', 'w'):
            raise ValueError(
                "Invalid mode string ''%s''. Only "
                "'a' and 'w' are acceptable modes " % mode)

        file = tables.open_file(filename, mode, title=title)

        # create the high-level structure of the cuds file
        for group in ('particle_container', 'lattice', 'mesh'):
            if "/" + group not in file:
                file.create_group('/', group, group)

        return cls(file)

    def close(self):
        """Closes a file

        """
        self._file.close()

    def add_lattice(self, name, lattice=None, keywords=[]):
        """Add lattice to the file.

        Parameters
        ----------
        name : str
            name of lattice
        lattice : Lattice, optional
            lattice to be added. If none is given,
            then an empty lattice is added.
        keywords : list of str, optional
            keywords to be copied from lattice to FileLattice.
            If 'none' then copy everything.


        Returns
        ----------
        FileLattice
            The lattice newly added to the file.  See
            get_lattice for more information.

        """
        if name in self._file.root.lattice:
            raise ValueError(
                'Lattice \'{n}\` already exists'.format(n=name))

        lat = FileLattice(self._file, name, lattice)
        self._lattices[name] = (lat, self._file.root.lattice)

        if lattice:
            # copy the contents of the lattice to the file
            for lattice_node in lattice.iter_nodes():
                # lat.add_node(lattice_node,keywords)
                lat.add_node(lattice_node, lattice_node.data.keys())

        self._file.flush()
        return lat

    def delete_lattice(self, name):
        """Delete lattice from file.

        Parameters
        ----------
        name : str
            name of lattice to delete
        """
        if name in self._lattices:
            self._file.getNode(self._lattices[name][1], name)._f_remove()
            del self._lattices[name]
        else:
            raise ValueError(
                'Lattice \'{n}\` does not exist'.format(n=name))

    def add_particle_container(self, name, particle_container=None):
        """Add particle container to the file.

        Parameters
        ----------
        name : str
            name of particle container
        particle_container : ABCParticleContainer, optional
            particle container to be added. If none is give,
            then an empty particle container is added.

        Returns
        ----------
        FileParticleContainer
            The particle container newly added to the file.  See
            get_particle_container for more information.

        """
        if name in self._file.root.particle_container:
            raise ValueError(
                'Particle container \'{n}\` already exists'.format(n=name))

        group = self._file.create_group('/particle_container/', name)
        pc = FileParticleContainer(group, self._file)
        self._particle_containers[name] = (pc, group)

        if particle_container:
            # copy the contents of the particle container to the file
            for particle in particle_container.iter_particles():
                pc.add_particle(particle)
            for bond in particle_container.iter_bonds():
                pc.add_bond(bond)

        self._file.flush()
        return pc

    def get_lattice(self, name):
        """Get lattice from file.

        The returned lattice can be used to query
        and change the related data stored in the file. If the
        file has been closed then the lattice should
        no longer be used.

        Parameters
        ----------
        name : str
            name of lattice to return
        """
        if name in self._lattices:
            return self._lattices[name][0]
        elif name in self._file.root.lattice:
            lat = FileLattice(self._file, name)
            self._lattices[name] = lat
            return lat
        else:
            raise ValueError(
                'Lattice \'{n}\` does not exist'.format(n=name))

    def get_particle_container(self, name):
        """Get particle container from file.

        The returned particle container can be used to query
        and change the related data stored in the file. If the
        file has been closed then the particle container should
        no longer be used.

        Parameters
        ----------
        name : str
            name of particle container to return
        """
        if name in self._particle_containers:
            return self._particle_containers[name][0]
        elif name in self._file.root.particle_container:
            group = tables.Group(self._file.root.particle_container, name)
            pc = FileParticleContainer(group, self._file)
            self._particle_containers[name] = pc
            return pc
        else:
            raise ValueError(
                'Particle container \'{n}\` does not exist'.format(n=name))

    def delete_particle_container(self, name):
        """Delete particle container from file.

        Parameters
        ----------
        name : str
            name of particle container to delete
        """
        if name in self._particle_containers:
            self._particle_containers[name][1]._f_remove(recursive=True)
            del self._particle_containers[name]
        else:
            raise ValueError(
                'Particle container \'{n}\` does not exist'.format(n=name))

    def iter_particle_containers(self, names=None):
        """Returns an iterator over a subset or all
        of the particle containers. The iterator iterator yields
        (name, particlecontainer) tuples for each particle container
        contained in the file.

        Parameters
        ----------
        names : list of str
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will
            be iterated over.

        """
        names = copy.deepcopy(names)
        if names is None:
            names = self._particle_containers.keys()
        for name in names:
            yield self.get_particle_container(name), name
