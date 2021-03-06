"""
    Testing module for lattice data classes.
"""
import unittest
import numpy as np
import numpy.testing as np_test
from math import sqrt
from simphony.core.cuba import CUBA
import simphony.cuds.lattice as la


class LatticeNodeTestCase(unittest.TestCase):
    """Test case for LatticeNode class."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_construct_lattice_node_default(self):
        """Creation of a lattice node."""
        node = la.LatticeNode((0, 0))

        self.assertIsInstance(node, la.LatticeNode,
                              "Error: not a LatticeNode!")
        self.assertEqual(node.id, (0, 0))

    def test_construct_lattice_node_copy(self):
        """Creation of a lattice node (copy constructor)."""
        node_org = la.LatticeNode((0, 1))
        node_org.data[CUBA.DENSITY] = 1.5
        node_org.data[CUBA.VELOCITY] = (0.2, -0.1)

        node_new = la.LatticeNode((2, 3), node_org.data)

        self.assertIsInstance(node_new, la.LatticeNode,
                              "Error: not a LatticeNode!")
        self.assertEqual(node_new.id, (2, 3))
        self.assertEqual(node_new.data[CUBA.DENSITY], 1.5)
        self.assertEqual(node_new.data[CUBA.VELOCITY], (0.2, -0.1))


class LatticeTestCase(unittest.TestCase):
    """Test case for Lattic class."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_construct_lattice_make(self):
        """Creation of lattices using the factory functions."""
        hexag_lat = la.make_hexagonal_lattice('Lattice1', 0.1, (11, 21))
        square_lat = la.make_square_lattice('Lattice2', 0.2, (12, 22))
        rectang_lat = la.make_rectangular_lattice('Lattice3',
                                                  (0.3, 0.35), (13, 23))
        cubic_lat = la.make_cubic_lattice('Lattice4',
                                          0.4, (14, 24), (4, 5, 6))
        orthop_lat = la.make_orthorombicp_lattice('Lattice5',
                                                  (0.5, 0.54, 0.58),
                                                  (15, 25), (7, 8, 9))

        self.assertIsInstance(hexag_lat, la.Lattice, "Error: not a Lattice!")
        self.assertIsInstance(square_lat, la.Lattice, "Error: not a Lattice!")
        self.assertIsInstance(rectang_lat, la.Lattice, "Error: not a Lattice!")
        self.assertIsInstance(cubic_lat, la.Lattice, "Error: not a Lattice!")
        self.assertIsInstance(orthop_lat, la.Lattice, "Error: not a Lattice!")

        self.assertEqual(hexag_lat.name, 'Lattice1')
        self.assertEqual(square_lat.name, 'Lattice2')
        self.assertEqual(rectang_lat.name, 'Lattice3')
        self.assertEqual(cubic_lat.name, 'Lattice4')
        self.assertEqual(orthop_lat.name, 'Lattice5')

        self.assertEqual(hexag_lat.type, 'Hexagonal')
        self.assertEqual(square_lat.type, 'Square')
        self.assertEqual(rectang_lat.type, 'Rectangular')
        self.assertEqual(cubic_lat.type, 'Cubic')
        self.assertEqual(orthop_lat.type, 'OrthorombicP')

        np_test.assert_array_equal(hexag_lat.size, (11, 21))
        np_test.assert_array_equal(square_lat.size, (12, 22))
        np_test.assert_array_equal(rectang_lat.size, (13, 23))
        np_test.assert_array_equal(cubic_lat.size, (14, 24))
        np_test.assert_array_equal(orthop_lat.size, (15, 25))

        np_test.assert_array_equal(hexag_lat.origin, (0, 0))
        np_test.assert_array_equal(square_lat.origin, (0, 0))
        np_test.assert_array_equal(rectang_lat.origin, (0, 0))
        np_test.assert_array_equal(cubic_lat.origin, (4, 5, 6))
        np_test.assert_array_equal(orthop_lat.origin, (7, 8, 9))

        np_test.assert_array_equal(hexag_lat.base_vect,
                                   (0.5*0.1, 0.5*sqrt(3)*0.1))
        np_test.assert_array_equal(square_lat.base_vect, (0.2, 0.2))
        np_test.assert_array_equal(rectang_lat.base_vect, (0.3, 0.35))
        np_test.assert_array_equal(cubic_lat.base_vect, (0.4, 0.4, 0.4))
        np_test.assert_array_equal(orthop_lat.base_vect, (0.5, 0.54, 0.58))

    def test_set_get_iter_lattice_nodes(self):
        """Creation of lattices using the factory functions."""
        rect_lat = la.make_rectangular_lattice('Lattice1', (0.1, 0.2),
                                               (10, 10), (1.5, 1.5))

        for i in range(10):
            node = rect_lat.get_node((i, i))
            node.data[CUBA.LABEL] = i
            rect_lat.update_node(node)

        node_count = 0
        for node in rect_lat.iter_nodes():
            if CUBA.LABEL in node.data:
                node_count += 1

        self.assertEqual(node_count, 10)

        iter = 0
        node_coords = np.zeros((10, 2), dtype=np.int32)
        for node in rect_lat.iter_nodes():
            if CUBA.LABEL in node.data:
                node_coords[iter, :] = node.id
                iter += 1

        iter = 0
        check_sum1 = 0
        for node in rect_lat.iter_nodes(node_coords):
            check_sum1 += node.data[CUBA.LABEL]
            coord = rect_lat.get_coordinate(node.id)

            self.assertEqual(node.id[0], iter)
            self.assertEqual(node.id[1], iter)
            self.assertEqual(coord[0], 0.1*iter + 1.5)
            self.assertEqual(coord[1], 0.2*iter + 1.5)

            iter += 1

        self.assertEqual(check_sum1, 45)

if __name__ == '__main__':
    unittest.main()
