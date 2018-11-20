import sys
import os
from collections import defaultdict
from flib.core.omim import OMIM


class Labels:

    def __init__(self, labels_dir=None, gmt=None, pos_label='1', neg_label='-1'):
        self._standards = {}

        if labels_dir:
            for labels_file in os.listdir(labels_dir):
                pos_genes, neg_genes = set(), set()
                with open(labels_dir + '/' + labels_file) as labelf:
                    lines = labelf.readlines()
                    for l in lines:
                        gene, label = l.strip('\t').split()[:2]
                        if label == pos_label:
                            pos_genes.add(gene)
                        elif label == neg_label:
                            neg_genes.add(gene)
                self._standards[labels_file] = (pos_genes, neg_genes)
        else:
            for gsid, genes in gmt.genesets.iteritems():
                pos_genes = genes
                neg_genes = gmt.genes - pos_genes
                self._standards[gsid] = (pos_genes, neg_genes)

    def get_labels(self, term_id):
        return self._standards[term_id]

    def get_terms(self):
        return self._standards.keys()


class OntoLabels:

    def __init__(self, obo=None, slim_terms=None):
        self._slim_terms = slim_terms
        self._obo = obo

    def get_labels(self, term_id):
        term = self._obo.get_term(term_id)
        if not term:
            return (set(), set())

        # Positive genes are annotated to the term
        pos = set(term.get_annotated_genes())

        unknown, all_genes = set(), set()
        term_tree = self._obo.get_ancestors(term_id)

        for obo_term in self._obo.get_termobject_list():
            obo_term_tree = self._obo.get_ancestors(obo_term.go_id)
            genes = set(obo_term.get_annotated_genes())

            # Terms share a slim in their ancestor tree
            if len(set(self._slim_terms & term_tree & obo_term_tree)):
                unknown |= genes
            all_genes |= genes

        neg = all_genes - unknown - pos

        return (pos, neg)


if __name__ == '__main__':
    from onto import DiseaseOntology

    do = DiseaseOntology.generate()
    OMIM().load_onto(onto=do)
    do.propagate()

    lines = open('../../files/do_slim.txt').readlines()
    slim_terms = set([l.strip() for l in lines])

    ol = OntoLabels(obo=do, slim_terms=slim_terms)
    (pos1, neg1) = ol.get_labels('DOID:0060041')
    print("%d, %d" %(len(pos1), len(neg1)))
