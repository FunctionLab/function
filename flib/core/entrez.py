import sys
import logging

import urllib2
import gzip
import io
from collections import defaultdict
from idmap import IDMap

logging.basicConfig()
logger = logging.getLogger(__name__)

URLS = {
    'homo_sapiens':
    'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
}

ENTREZID, SYMBOL, XREFS = 1, 2, 5


class Entrez:

    def __init__(self):
        self._symbols = defaultdict(set)
        self._xrefs = defaultdict(set)
        self._entrez = defaultdict(set)

    def load(self, organism='homo_sapiens'):

        gene_zip_fh = urllib2.urlopen(
            URLS[organism], timeout=5)
        geneinfo = gzip.GzipFile(fileobj=io.BytesIO(gene_zip_fh.read()))
        gene_zip_fh.close()

        for line in geneinfo:
            toks = line.strip().split('\t')
            symbol, entrezid, xrefs = toks[SYMBOL], toks[ENTREZID], toks[XREFS]

            self._symbols[symbol].add(entrezid)
            self._entrez[entrezid].add(symbol)

            for xref in xrefs.split('|'):
                self._xrefs[xref].add(entrezid)
                self._entrez[entrezid].add(xref)

    def get(self, symbol=None, xref=None):
        if symbol and symbol in self._symbols:
            return self._symbols[symbol]
        elif xref and xref in self._xrefs:
            return self._xrefs[xref]
        else:
            return None

    def get_symbol_map(self):
        return IDMap(key_map = self._symbols)

    def get_xref_map(self):
        return IDMap(key_map = self._xrefs)

if __name__ == '__main__':
    entrez = Entrez()
    entrez.load()
    print entrez.get(symbol='LRP6')
