import logging
import requests
import urllib2

from onto import GeneOntology
from entrez import Entrez
from flib import settings

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class GOA:
    '''Gene ontology associations'''

    def __init__(self, org='Homo sapiens'):
        self._onto = None
        self._org = org
        self._meta = {}

    def load_onto(self, onto=None, idmap=None):
        if not onto:
            onto = GeneOntology.generate()

        for prefix, suffix in zip(settings.GOA_PREFIX, settings.GOA_ASSOC_SUFFIX):
            annot_zip = settings.GOA_ASSOC_URL + \
                ''.join((prefix, settings.GOA_NAMES[self._org], suffix))
            ret = requests.head(annot_zip)
            if ret.status_code < 400:
                logger.info('Loading: %s', annot_zip)
                onto.populate_annotations(
                    annot_zip,
                    remote_location=True)
                break
            else:
                logger.debug('URL not available: %s', annot_zip)

        for prefix, suffix in zip(settings.GOA_PREFIX, settings.GOA_INFO_SUFFIX):
            info = settings.GOA_ASSOC_URL + \
                ''.join((prefix, settings.GOA_NAMES[self._org], suffix))
            ret = requests.head(info)
            if ret.status_code < 400:
                logger.info('Loading: %s', info)
                annot_info = urllib2.urlopen(info, timeout=5)
                annot_info = eval(annot_info.read())
                self._meta = annot_info
                break
            else:
                logger.debug('URL not available: %s', annot_zip)

        if idmap:
            onto.map_genes(idmap, xdb_prefixed=True)

        self._onto = onto
        return onto

    def get_meta_data(self, key):
        return self._meta.get(key)


if __name__ == '__main__':
    entrez_map = Entrez()
    entrez_map.load()

    goa = GOA()
    onto = goa.load_onto()
    onto.print_to_gmt_file('go.gmt')
