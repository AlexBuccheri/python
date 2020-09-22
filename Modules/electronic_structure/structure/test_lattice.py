import unittest
import numpy as np

from modules.electronic_structure.structure import lattice
from modules.electronic_structure.structure.bravais import simple_cubic, body_centred_cubic, face_centred_cubic

class MyTestCase(unittest.TestCase):
    """ Unit tests for lattice.py module """


    def test_parallelpiped_volume_simple_cubic(self):
        """ Test volume for simple cubic lattice"""
        al = 5
        lattice_vectors = simple_cubic(5)
        vol = lattice.parallelpiped_volume(lattice_vectors)
        self.assertEqual(vol, al**3)


    def test_reciprocal_lattice_vectors(self):
        """
        a_i . b_j = 2pi * delta_ij
        where {a} are real-space lattice vectors and
        {b} are reciprocal-space lattice vectors
        """
        from modules.maths.matrix_utils import off_diagonal_indices

        a = face_centred_cubic(5.1)
        b = lattice.reciprocal_lattice_vectors(a)
        self.assertTrue(a.shape == (3,3))
        self.assertTrue(b.shape == (3, 3))
        dp = np.empty(shape=(3, 3))

        for i in range(0, 3):
            for j in range(0, 3):
                dp[i, j] = np.dot(a[:, i], b[:, j])

        off_diagonal_index = off_diagonal_indices(np.size(dp, 0))
        off_diagonals = dp[off_diagonal_index]

        self.assertTrue(np.all(dp.diagonal() == 2 * np.pi))
        self.assertTrue(np.all(off_diagonals == 0))


    def test_translation_integers_for_simple_cubic_cell(self):
        """
        For a cubic system, expect the naive expression
        ceil(cutoff / |a|) to give the same maximum integer as
        translation_integers_for_radial_cutoff function,
        where a is the cubic lattice vector
        """
        lattice_vectors = simple_cubic(5)
        cutoff = 21
        integers = lattice.simple_cubic_cell_translation_integers(lattice_vectors, cutoff)
        integers2 = lattice.translation_integers_for_radial_cutoff(lattice_vectors, cutoff)
        self.assertTrue(np.all(integers == integers2))


    def get_extremal_translations(self, lattice_vectors, max_integers):
        """ Used by tests but not a test itself"""
        extremal_integers = lattice.get_extremal_integers(max_integers[0],
                                                          max_integers[1],
                                                          max_integers[2])
        self.assertTrue(extremal_integers.shape == (3,26))

        extremal_translations = np.matmul(lattice_vectors, extremal_integers)
        extremal_translation_magnitudes = np.empty(shape=(26))
        for i in range(0, 26):
            extremal_translation_magnitudes[i] = \
                np.linalg.norm(extremal_translations[:, i])

        return extremal_translation_magnitudes


    def test_translation_integers_from_radial_cutoff(self):
        """
        For non-cubic cells, make sure the extremial
        translation vectors |{T}_ext|, defined by all permutations of
        the max lattice integers (n, m, l), are all less than or equal
        to the radial translation cutoff
        """

        def test_simple_cubic(cutoff):
            lattice_vectors = simple_cubic(5)

            #TODO ALEX FIX THIS CALL
            max_integers = \
                lattice.translation_integers_for_radial_cutoff(lattice_vectors, cutoff)
            self.assertEqual(max_integers[0], 5)
            self.assertEqual(max_integers[1], 5)
            self.assertEqual(max_integers[2], 5)

            extremal_translation_magnitudes = \
                self.get_extremal_translations(lattice_vectors, max_integers)

            # For all extremal translations, is |T| < cutoff ?
            # I.e. is the set of integers not large enough for the radial cutoff?
            indices = np.where(extremal_translation_magnitudes < cutoff)[0]
            self.assertTrue(indices.size == 0)

        def test_bcc(cutoff):
            lattice_vectors = body_centred_cubic(5)
            #TODO ALEX FIX THIS CALL
            max_integers = lattice.translation_integers_for_radial_cutoff(lattice_vectors,
                                                                          cutoff)
            self.assertEqual(max_integers[0], 6)
            self.assertEqual(max_integers[1], 6)
            self.assertEqual(max_integers[2], 6)

            extremal_translation_magnitudes = self.get_extremal_translations(lattice_vectors,
                                                                             max_integers)
            indices = np.where(extremal_translation_magnitudes < cutoff)[0]
            self.assertTrue(indices.size == 0)

        def test_fcc(cutoff):
            lattice_vectors = face_centred_cubic(5)

            naive_max_integers = lattice.simple_cubic_cell_translation_integers(lattice_vectors,
                                                                                cutoff)

            max_integers = lattice.another_translation_integers_for_radial_cutoff(lattice_vectors,
                                                                          cutoff)
            self.assertEqual(max_integers[0], 8)
            self.assertEqual(max_integers[1], 8)
            self.assertEqual(max_integers[2], 8)
            self.assertTrue(np.all(naive_max_integers != max_integers))

            def lattice_sum(lattice, n, cutoff):
                cutoff_squared = cutoff * cutoff
                translations = []
                for k in range(-n[2], n[2]):
                    for j in range(-n[1], n[1]):
                        for i in range(-n[0], n[0]):
                            r = np.matmul(lattice, np.array([i, j, k]))
                            if np.dot(r, r) <= cutoff_squared:
                                translations.append(list(r))
                return translations

            less_translations = lattice_sum(lattice_vectors, naive_max_integers, cutoff)
            translations = lattice_sum(lattice_vectors, max_integers, cutoff)
            print(set(translations)-set(less_translations))

            naive_extremal_translation_magnitudes = \
                self.get_extremal_translations(lattice_vectors, naive_max_integers)
            extremal_translation_magnitudes = \
                self.get_extremal_translations(lattice_vectors, max_integers)

            naive_indices = np.where(cutoff > naive_extremal_translation_magnitudes)[0]
            indices = np.where(cutoff > extremal_translation_magnitudes)[0]

            # This first one is surprising. I would expect a few extremal translation magnitudes
            # to be less than the cutoff, based on my visual assertment of the supercell (i.e. it doesn't look spherical)
            # MUST be the smallest ones that are ultimately the problem
            print(naive_extremal_translation_magnitudes)

            # NO IDEA why this isn't working. Could plot out all the points on the surface of a sphere

            # TODO(Alex). Could sample the surface of a surface, and the surface of this cell,
            # and check that all points on the sphere are within the cell - seems like overkill
            self.assertTrue(naive_indices.size == 0)
            self.assertTrue(indices.size == 0)


        # Radial cut-off defined w.r.t. central cell
        cutoff = 21
        test_simple_cubic(cutoff)
        test_bcc(cutoff)
        test_fcc(cutoff)
        # Should in principle test for all bravais lattices but
        # intend to do this in entos instead




if __name__ == '__main__':
    unittest.main()
