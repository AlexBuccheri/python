import unittest
import numpy as np

import cell_operations

class MyTestCase(unittest.TestCase):

    directories = cell_operations.Directories(structure='aei',
                                              input='inputs',
                                              output='aei_outputs')

    unit_cell, lattice_vectors = cell_operations.get_primitive_unit_cell(
        directories, visualise=True)

    translations = cell_operations.translations_for_fully_coordinated_unit(
        unit_cell, lattice_vectors)

    def test_get_primitive_unit_cell(self):
        lattice = np.array([[18.5277,      0.,      0.],
                            [     0.,  6.8556, -6.8556],
                            [     0.,  6.3323,  6.3323]])
        self.assertEqual(len(self.unit_cell), 72)
        self.assertTrue(np.array_equal(self.lattice_vectors, lattice))

        species = [atom.species for atom in self.unit_cell]
        self.assertTrue(all(element == 'Si' or element == 'O'
                            for element in species))

    def test_translations_for_fully_coordinated_unit(self):
        self.assertEqual(len(self.translations), 27)
        middle_index = int(0.5 * len(self.translations))
        translation = self.translations[middle_index].tolist()
        self.assertListEqual(translation, [0, 0, 0],
                             "Expect central translation to be [0,0,0]")

    # #This currently doesn't work - might be same atoms in different order
    # def test_replace_uncoordinated_atoms(self):
    #     """
    #     Test folding replacement atoms back into the central cell
    #     Should retrieve the original cell
    #     """
    #     # Unit cell with some atoms outside of the central unit
    #     unit_cell_replaced = cell_operations.replace_uncoordinated_atoms(
    #         self.unit_cell, self.translations)
    #     unit_cell_restored = cell_operations.ensure_atoms_in_central_cell(
    #         unit_cell_replaced, self.lattice_vectors)
    #
    #     # Can assert per atom
    #     for ia in range(0, len(self.unit_cell)):
    #         print(self.unit_cell[ia].position, unit_cell_restored[ia].position)




if __name__ == '__main__':
    unittest.main()
