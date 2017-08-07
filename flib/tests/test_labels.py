import unittest
import numpy

from flib.core.obo import OBO, GOTerm
from flib.core.onto import GeneOntology
from flib.core.labels import OntoLabels

GO_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/gene-associations/gene_association.goa_human.gz?rev=1.353'
DSREPAIR_ID = 'GO:0000724'
DNA_REP_ID = 'GO:0006260'
DNA_METAB = 'GO:0006259'

class TestLabels(unittest.TestCase):

    def setUp(self):
        self.go = OBO('files/test_data/go.obo')
        self.go.populate_annotations(GO_URL, remote_location=True)

        self.dsrepair_term = self.go.get_term(DSREPAIR_ID)

        lines = open('files/go_neg_slim.txt').readlines()
        self.slim_terms = set([l.strip() for l in lines])

    def testOntoLabelsNeg(self):
        ol = OntoLabels(obo=self.go, slim_terms=self.slim_terms)

        (pos, neg) = ol.get_labels(DSREPAIR_ID)
        self.assertTrue(len(pos) > 0)
        self.assertEqual(pos,set(self.dsrepair_term.get_annotated_genes()))

        similar_term = self.go.get_term(DNA_REP_ID).get_annotated_genes()

        slim_overlap = self.go.get_ancestors(DNA_REP_ID) & \
            self.go.get_ancestors(DSREPAIR_ID) & \
            self.slim_terms

        self.assertTrue(DNA_METAB in slim_overlap)
        self.assertTrue(len(neg & set(similar_term)) == 0)

        for dterm in self.go.get_descendents(DNA_METAB):
            dgenes = set(self.go.get_term(dterm).get_annotated_genes())
            self.assertTrue(len(neg & dgenes) == 0)

