import logging
import requests
#import urllib2

from flib.core.onto import GeneOntology
from flib.core.url import URLResource
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

            logger.info("Loading annotations from: "+ annot_zip)
            onto.populate_annotations(annot_zip, remote_location=True)


        for prefix, suffix in zip(settings.GOA_PREFIX, settings.GOA_INFO_SUFFIX):
            info = settings.GOA_ASSOC_URL + \
                ''.join((prefix, settings.GOA_NAMES[self._org], suffix))

            logger.info("Loading from: "+ info)
            annot_info = URLResource(info).get_lines()
            self._meta = annot_info


        if idmap:
            onto.map_genes(idmap, xdb_prefixed=True)

        self._onto = onto
        return onto

    def get_meta_data(self, key):
        return self._meta.get(key)


if __name__ == '__main__':
    # entrez_map = Entrez()
    # entrez_map.load()

    goa = GOA()
    onto = goa.load_onto()
    print(onto)
    onto.print_to_gmt_file('go.gmt')
