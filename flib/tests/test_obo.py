import unittest
import numpy

from flib.core.obo import OBO, GOTerm

GO_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/ontology/gene_ontology.obo?rev=4.2451'
DSREPAIR_ID = 'GO:0000724'

class TestOBO(unittest.TestCase):

    def setUp(self):
        self.go = OBO('files/test_data/go.obo')
        self.dsrepair_term = self.go.get_term(DSREPAIR_ID)

    def test_load_remote(self):
        '''Test loading an obo file from a URL'''
        self.go = OBO()
        self.go.load_obo(GO_URL, remote_location=True)
        self.test_load()

    def test_load(self):
        '''Test an obo is loaded'''
        # Check that root biologcal process term is loaded
        term = self.go.get_term('GO:0008150')
        self.assertTrue(term is not None and term.name == 'biological_process')
        self.assertTrue(len(term.parent_of) > 0)

    def test_add_annotation(self):
        '''Test that adding a gene annotations is correct'''
        self.assertTrue(self.dsrepair_term is not None)

        self.dsrepair_term.add_annotation('672')
        genes = self.dsrepair_term.get_annotated_genes()
        self.assertEqual(genes[0],'672')

    def test_propagation(self):
        '''Test that gene propagation adds annotations to the correct terms'''
        self.assertTrue(self.dsrepair_term is not None)

        self.dsrepair_term.add_annotation('672')
        self.go.propagate()
        term_count = 0
        for term in self.go.get_termobject_list():
            term_count += len(term.get_annotated_genes())

        self.assertEqual(term_count, 26)

    def test_ancestors(self):
        '''Test that ancestor terms are correct'''
        self.assertTrue(self.dsrepair_term is not None)
        self.assertEqual(len(self.go.get_ancestors(self.dsrepair_term.go_id)), 25)

    def test_parents(self):
        '''Test that parent terms are stored correctly'''
        parents = self.dsrepair_term.child_of
        self.assertEqual(len(parents),2)
        self.assertTrue(GOTerm('GO:0006302') in parents)
        self.assertTrue(GOTerm('GO:0000725') in parents)

    def test_term_equals(self):
        '''Test GOTerm equals method'''
        self.assertEqual(self.dsrepair_term, GOTerm(DSREPAIR_ID) )

    def tearDown(self):
        self.go = None
        self.dsrepair_term = None
