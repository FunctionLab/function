import unittest
import numpy

from flib.core.obo import OBO

GO_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/ontology/gene_ontology.obo?rev=4.2451'

class TestOBO(unittest.TestCase):

    def setUp(self):
        self.go = OBO('files/test_data/go.obo')

    def test_load_remote(self):
        self.go = OBO()
        self.go.load_obo(GO_URL, remote_location=True)
        self.test_load()

    def test_load(self):
        # Check that root biologcal process term is loaded
        term = self.go.get_term('GO:0008150')
        self.assertTrue(term is not None and term.name == 'biological_process')
        self.assertTrue(len(term.parent_of) > 0)

    def test_add_annotation(self):
        term = self.go.get_term('GO:0000724')
        self.assertTrue(term is not None)

        term.add_annotation('672')
        genes = term.get_annotated_genes()
        self.assertEqual(genes[0],'672')

    def test_propagation(self):
        dsrepair_term = self.go.get_term('GO:0000724')
        self.assertTrue(dsrepair_term is not None)

        dsrepair_term.add_annotation('672')
        self.go.propagate()
        term_count = 0
        for term in self.go.get_termobject_list():
            term_count += len(term.get_annotated_genes())

        self.assertEqual(term_count, 26)

    def test_ancestors(self):
        dsrepair_term = self.go.get_term('GO:0000724')
        self.assertTrue(dsrepair_term is not None)
        self.assertEqual(len(self.go.get_ancestors(dsrepair_term.go_id)), 25)

    def tearDown(self):
        self.go = None
