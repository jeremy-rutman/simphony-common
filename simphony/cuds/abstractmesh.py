"""
Module for Abstract Mesh class

"""

from abc import ABCMeta, abstractmethod


class ABCMesh(object):
    """Abstract base class for mesh."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def get_point(self, uuid):
        pass

    @abstractmethod
    def get_edge(self, uuid):
        pass

    @abstractmethod
    def get_face(self, uuid):
        pass

    @abstractmethod
    def get_cell(self, uuid):
        pass

    @abstractmethod
    def add_point(self, point):
        pass

    @abstractmethod
    def add_edge(self, edge):
        pass

    @abstractmethod
    def add_face(self, face):
        pass

    @abstractmethod
    def add_cell(self, cell):
        pass

    @abstractmethod
    def update_point(self, point):
        pass

    @abstractmethod
    def update_edge(self, edge):
        pass

    @abstractmethod
    def update_face(self, face):
        pass

    @abstractmethod
    def update_cell(self, cell):
        pass

    @abstractmethod
    def iter_points(self, point_ids=[]):
        pass

    @abstractmethod
    def iter_edges(self, edge_ids=[]):
        pass

    @abstractmethod
    def iter_faces(self, face_ids=[]):
        pass

    @abstractmethod
    def iter_cells(self, cell_ids=[]):
        pass

    @abstractmethod
    def has_edges(self):
        pass

    @abstractmethod
    def has_faces(self):
        pass

    @abstractmethod
    def has_cells(self):
        pass
