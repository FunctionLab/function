import unittest
import numpy

from flib.core.obo import OBO

GO_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/ontology/gene_ontology.obo?rev=4.2451'


class TestOBO(unittest.TestCase):

    def setUp(self):
        self.go = OBO()

    def test_load_remote(self):
        self.go.load_obo(GO_URL, remote_location=True)

        # Check that root biologcal process term is loaded
        term = self.go.get_term('GO:0008150')
        self.assertTrue(term is not None and term.name == 'biological_process')
        self.assertTrue(len(term.parent_of) > 0)

    def tearDown(self):
        self.go = None
