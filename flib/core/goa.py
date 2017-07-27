import sys
import logging
from collections import defaultdict
import re
import requests
import urllib2
import gzip
import io

from onto import GeneOntology
from entrez import Entrez
from idmap import IDMap

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

GO_NAMES = {
    'Arabidopsis thaliana': 'tair',
    'Homo sapiens': 'human',
    'Mus musculus': 'mgi',
    'Rattus norvegicus': 'rgd',
    'Danio rerio': 'zfin',
    'Drosophila melanogaster': 'fb',
    'Saccharomyces cerevisiae': 'sgd',
    'Caenorhabditis elegans': 'wb',
    'Pseudomonas aeruginosa': 'pseudocap'
}

GO_PREFIX = ['goa_', 'gene_association.']
GO_ASSOC_SUFFIX = ['.gaf.gz', '.gz']
GO_INFO_SUFFIX = ['.gaf.json', '.json']

GO_ASSOC_URL = 'http://www.geneontology.org/gene-associations/'
GO_VERSION_KEY = 'submissionDate'

class GOA:
    '''Gene ontology associations'''
    def __init__(self, org = 'Homo sapiens'):
        self._onto = None
        self._org = org

    def load_onto(self, onto=None, idmap=None):
        if not onto:
            onto = GeneOntology.generate()

        for prefix, suffix in zip(GO_PREFIX, GO_ASSOC_SUFFIX):
            annot_zip = GO_ASSOC_URL + \
                ''.join((prefix, GO_NAMES[self._org], suffix))
            ret = requests.head(annot_zip)
            if ret.status_code < 400:
                logger.info('Loading: %s', annot_zip)
                onto.populate_annotations(
                    annot_zip,
                    remote_location=True)
                break
            else:
                logger.info('URL not available: %s', annot_zip)

        if idmap:
            onto.map_genes(idmap, xdb_prefixed=True)

        self._onto = onto
        return onto

if __name__ == '__main__':
    entrez_map = Entrez()
    entrez_map.load()

    goa = GOA()
    onto = goa.load_onto(idmap=entrez_map.get_xref_map())
    onto.print_to_gmt_file('go.gmt')
