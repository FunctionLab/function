import logging
from collections import defaultdict
from six import iteritems

from flib.core.onto import DiseaseOntology
from flib.core.entrez import Entrez
from flib.core.url import URLResource
from flib import settings
logging.basicConfig()
logger = logging.getLogger(__name__)

DISEASE, REPORTED_GENES, MAPPED_GENES, MAPPED_TRAIT, MAPPED_TRAIT_URI = 7, 13, 14, 34, 35
COLS = [DISEASE, REPORTED_GENES, MAPPED_GENES, MAPPED_TRAIT, MAPPED_TRAIT_URI]

class GWASCatalog:

    def __init__(self):
        self._onto = None
        self._data = None

    def load_onto(self, onto=None, idmap=None):
        if not self._data:
            self.load_data()

        if not onto:
            onto = DiseaseOntology.generate()

        xrefs = onto.get_xref_mapping('EFO')

        for (trait, uri), genes in iteritems(self._data):
            uids = [x.split('/')[-1].replace('EFO_', '')
                    for x in uri.split(',')]
            if len(uids) > 1:
                logger.info('Multiple mappings %s', uri)
                continue

            terms = [onto.get_term(termid)
                     for uid in uids for termid in xrefs[uid]]
            for term in terms:
                for gene in genes:
                    mapped_genes = idmap[gene] if idmap else (gene,)
                    for gid in mapped_genes:
                        term.add_annotation(gid=gid)

        self._onto = onto
        return onto

    def load_data(self):
        lines = URLResource(settings.GWAS_URL).get_lines()

        genesets = defaultdict(set)
        headers = []
        for (i, line) in enumerate(lines):
            tok = line.strip().split('\t')
            if not i:
                headers = tok
            else:
                if len(tok) <= COLS[-1]:
                    logger.error('Error on line %i', i)
                    continue

                key = (tok[MAPPED_TRAIT], tok[MAPPED_TRAIT_URI])
                genes = set([x.strip()
                             for x in tok[REPORTED_GENES].split(',')])
                genesets[key] |= genes

        self._data = genesets

        return True

if __name__ == '__main__':
    entrez_map = Entrez()
    entrez_map.load()

    gwas = GWASCatalog()
    onto = gwas.load_onto(idmap=entrez_map.get_symbol_map())

    onto.print_to_gmt_file('test.txt')
