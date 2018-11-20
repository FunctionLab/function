from __future__ import absolute_import
import sys
import logging

from collections import defaultdict

from flib.core.idmap import IDMap
#from .idmap import IDMap
from flib import settings
from flib.core.url import URLResource

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

ENTREZID, SYMBOL, XREFS = 1, 2, 5

class Entrez:

    def __init__(self):
        self._symbols = defaultdict(set)
        self._xrefs = defaultdict(set)
        self._entrez = defaultdict(set)

    def load(self, organism='Homo sapiens'):

        logger.info("Loading from "+settings.GENEINFO_URLS[organism])
        geneinfo = URLResource(settings.GENEINFO_URLS[organism]).get_lines()

        for line in geneinfo:
            toks = line.strip().split('\t')
            if len(toks) == 1:
                continue
            symbol, entrezid, xrefs = toks[SYMBOL], toks[ENTREZID], toks[XREFS]

            self._symbols[symbol].add(entrezid)
            self._entrez[entrezid].add(symbol)

            for xref in xrefs.split('|'):
                self._xrefs[xref].add(entrezid)
                self._entrez[entrezid].add(xref)

        # Add uniprot IDs
        uniprot = URLResource(settings.UNIPROT_URLS[organism]).get_lines()

        for line in uniprot:
            tok = line.strip().split('\t')
            if len(tok) == 1:
                continue
            if len(tok) == 3 and tok[1] == 'GeneID':
                uniprotid, itype, entrezid = tok
                uniprotid = settings.UNIPROT_PREFIX + ':' + uniprotid
                self._xrefs[uniprotid].add(entrezid)
                self._entrez[entrezid].add(uniprotid)


    def get(self, symbol=None, xref=None):
        if symbol and symbol in self._symbols:
            return self._symbols[symbol]
        elif xref and xref in self._xrefs:
            return self._xrefs[xref]
        else:
            return None

    def get_symbol_map(self):
        return IDMap(key_map=self._symbols)

    def get_xref_map(self):
        return IDMap(key_map=self._xrefs)

    def get_entrez_map(self):
        return IDMap(key_map=self._entrez)

if __name__ == '__main__':
    entrez = Entrez()
    entrez.load()
    #print(entrez.get_symbol_map().keys())
