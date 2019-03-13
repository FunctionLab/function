import logging

from flib.core.onto import GeneOntology
from flib.core.url import URLResource
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

            logger.info("Loading meta info from: "+ info)
            annot_info = URLResource(info).get_lines()
            print(annot_info)
            self._meta = annot_info

        if idmap:
            onto.map_genes(idmap, xdb_prefixed=True)

        self._onto = onto
        return onto

    def get_meta_data(self, key):
        print((self._meta))
        return self._meta.get(key)


if __name__ == '__main__':
    # entrez_map = Entrez()
    # entrez_map.load()

    goa = GOA()
    onto = goa.load_onto()
    print((goa._meta))
    #onto.print_to_gmt_file('go.gmt')
    #onto.print_to_dir('test')
    #onto.print_to_single_file('go.single')
    #onto.print_to_mat_file('go.mat')
