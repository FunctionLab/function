import unittest
import numpy

from flib.core.dab import Dab

class TestDab(unittest.TestCase):

    def setUp(self):
        self.dab_file = 'files/test_data/test_dab.dab'
        self.dab = Dab(self.dab_file)

        self.dat_file = 'files/test_data/test_dat.dat'
        self.dat = open(self.dat_file).readlines()

    def tearDown(self):
        return

    def testOpenDab(self):
        # Test the total number of loaded genes
        self.assertEqual(len(self.dab.gene_list), 16)

        # Test the total number of values (16 choose 2)
        self.assertEqual(len(self.dab.dat), 120)

    def testGetValue(self):
        # Test the values from Dab class with values from dab exported as dat
        for l in self.dat:
            g1, g2, value = l.strip().split('\t')
            val = float(value)
            dat_val = self.dab.get_value_genestr(g1,g2)
            assert numpy.isclose(val, dat_val, rtol=1e-05, atol=1e-08)

    def testGet(self):
        for i,g1 in enumerate(self.dab.gene_list):
            vals = self.dab.get(g1)
            for j,g2 in enumerate(self.dab.gene_list):
                if g1 == g2:
                    # Self interaction should be 1
                    self.assertEqual(vals[j],1)
                else:
                    # Test the values from dab.get match dab.get_value
                    self.assertEqual(vals[j], self.dab.get_value_genestr(g1,g2))

