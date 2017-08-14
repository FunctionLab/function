import sys
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

import re
import requests

from flib.settings import OMIM_MIM2GENE, OMIM_GENEMAP, \
    OMIM_LIMIT_TYPE, OMIM_LIMIT_PHENO, OMIM_LIMIT_STATUS
from flib.core.onto import DiseaseOntology

# This should be standardized, but you never know
# Most disorders that have omimphenotypes fit this expression
FIND_MIMID = re.compile('\, [0-9]* \([1-4]\)')


class mim_disease:

    def __init__(self):
        self.mimid = ''
        self.is_susceptibility = 0  # Whether it has {}
        self.phe_mm = ''  # Phenotype mapping method
        self.genetuples = []  # (Gene ID, Gene Status)


class OMIM:

    def __init__(self, key=None):
        self._key = key
        self._onto = None
        self._data = None
        self._meta = {}

    def load_data(self):

        mim_gene = {}
        mim2gene_list = requests.get(OMIM_MIM2GENE).text.splitlines()

        for line in mim2gene_list:  # Loop from Dima @ Princeton
            if line.startswith('#'):
                continue
            toks = line.split('\t')
            if len(toks) < 3:
                logger.error('Can\'t parse line: %s', line)
                continue
            mim, gtype, gid = toks[:3]
            if gtype in OMIM_LIMIT_TYPE:
                if mim in mim_gene:
                    logger.warning("MIM already exists: %s", mim)
                if gid:
                    mim_gene[mim] = gid

        mimdiseases = {}
        genemap_list = requests.get(OMIM_GENEMAP).text.splitlines()
        genemap_version = None

        # TODO: Add support for publications
        for l in genemap_list:  # Loop from Dima @ Princeton
            if l.startswith('#'):
                if l.startswith('# Generated:'):
                    genemap_version = l.split('# Generated:')[1].strip()
                    self._meta['genemap_version'] = genemap_version
                continue

            l_split = l.split('\t')
            status = l_split[6].strip()
            mim_geneid = l_split[8].strip()
            disorders = l_split[11].strip()

            if disorders != '' and status in OMIM_LIMIT_STATUS and mim_geneid in mim_gene:
                logger.debug('%s with status %s', disorders, status)

                geneid = mim_gene[mim_geneid]
                tuple_gid_status = (geneid, status)

                disorders_list = disorders.split(';')
                for d in disorders_list:
                    if '[' not in d and '?' not in d:
                        mim_info = re.search(FIND_MIMID, d)
                        if mim_info:
                            # print 'Has necessary info'
                            # TODO: Make sure to include ? and [
                            info_split = mim_info.group(0).split(' ')
                            mim_disease_id = info_split[1].strip()
                            mim_phetype = info_split[2].strip()
                            if mim_phetype == OMIM_LIMIT_PHENO:
                                # print 'Correct phenotype'
                                if mim_disease_id not in mimdiseases:
                                    mimdiseases[mim_disease_id] = mim_disease()
                                    mimdiseases[
                                        mim_disease_id].mimid = mim_disease_id
                                    mimdiseases[
                                        mim_disease_id].phe_mm = mim_phetype
                                if '{' in d:
                                    mimdiseases[
                                        mim_disease_id].is_susceptibility = 1
                                if tuple_gid_status not in mimdiseases[
                                        mim_disease_id].genetuples:
                                    mimdiseases[mim_disease_id].genetuples.append(
                                        tuple_gid_status)

        self._data = mimdiseases
        return True

    def load_onto(self, onto=None, idmap=None):
        if not self._data:
            self.load_data()

        if not onto:
            onto = DiseaseOntology.generate()

        xrefs = onto.get_xref_mapping('OMIM')

        for (omim_id, mim_entry) in self._data.iteritems():
            if omim_id not in xrefs:
                continue

            d_or_s = 'S' if mim_entry.is_susceptibility else 'D'

            for doid in xrefs[omim_id]:
                term = onto.get_term(doid)
                for g in mim_entry.genetuples:

                    genes = idmap[g[0]] if idmap else (g[0],)

                    for gid in genes:
                        term.add_annotation(gid=gid)

        self._onto = onto
        return onto

    def get_meta_data(self, key):
        return self._meta.get(key)

if __name__ == '__main__':
    omim = OMIM()
    onto = omim.load_onto()

    #onto.print_to_gmt_file('test.txt')
    print omim._meta
