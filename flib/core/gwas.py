import sys
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
from collections import defaultdict

import re

import requests
import requests_ftp

from flib.core.do import DiseaseOntology

requests_ftp.monkeypatch_session()

GWAS_URL = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'
DISEASE, REPORTED_GENES, MAPPED_GENES, MAPPED_TRAIT, MAPPED_TRAIT_URI  = 7, 13, 14, 34, 35
COLS = [DISEASE, REPORTED_GENES, MAPPED_GENES, MAPPED_TRAIT, MAPPED_TRAIT_URI]

class GWASCatalog:

    def __init__(self):
        self._onto = None

    def load(self, onto=None):
        if not onto:
            onto = DiseaseOntology.generate()

        term_map = defaultdict(set)
        for (term_id, term) in onto.go_terms.iteritems():
            if 'EFO' in term.xrefs:
                for xref in term.xrefs['EFO']:
                    term_map['EFO_'+xref].add(term)

        lines = requests.get(GWAS_URL).text.encode('utf-8').splitlines()

        genesets, disease = defaultdict(set), defaultdict(set)
        headers = []
        for (i, line) in enumerate(lines):
            tok = line.strip().split('\t')
            if not i:
                headers = tok
            else:
                if len(tok) <= COLS[-1]:
                    logger.error('Error on line %i', i)
                    continue

                key = (tok[MAPPED_TRAIT],tok[MAPPED_TRAIT_URI])  
                genes = set([x.strip() for x in tok[REPORTED_GENES].split(',')])
                genesets[key] |= genes 
                disease[key].add(tok[DISEASE])        

        for (trait, uri), genes in genesets.iteritems():
            uids = [x.split('/')[-1] for x in uri.split(',')]
            terms = [term for uid in uids for term in term_map[uid]]
            if len(terms) != len(uids):
                logger.info('Multiple mappings %s', uri)
                continue

            for term in terms:
                for g in genes:
                    term.add_annotation(gid=g)

        self._onto = onto
        return onto

if __name__ == '__main__':
    gwas = GWASCatalog()
    onto = gwas.load()

    onto.print_to_gmt_file('test.txt')
