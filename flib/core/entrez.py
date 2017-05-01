import sys
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

import urllib2
import gzip
import io

import idmap



URLS = {
    'homo_sapiens' : 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'        
}

SYMBOL, ENTREZID, XREFS = 1, 2, 5

class Entrez:

    def __init__(self):
        return

    def load(self, organism='homo_sapiens'):

        gene_zip_fh = urllib2.urlopen(
            URLS[organism], timeout=5)
        geneinfo = gzip.GzipFile(fileobj=io.BytesIO(gene_zip_fh.read()))
        gene_zip_fh.close()
        
        for line in geneinfo:
            toks = line.strip().split('\t') 
            print toks[SYMBOL], toks[ENTREZID], toks[XREFS]


if __name__ == '__main__':
    entrez = Entrez()
    entrez.load()
