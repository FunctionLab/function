import sys
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

import re

import requests
import requests_ftp

from flib.core.do import DiseaseOntology

requests_ftp.monkeypatch_session()

MIM2GENE = 'http://omim.org/static/omim/data/mim2gene.txt'
GENEMAP = 'http://data.omim.org/downloads/NA8NpTI7QLK_CpW0PqV5uw/genemap.txt'

LIMIT_TYPE = set(['gene', 'gene/phenotype'])
LIMIT_PHENO = '(3)'
LIMIT_STATUS = ['C', 'P']

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

    def load(self, onto=None):
        if not onto:
            onto = DiseaseOntology.generate()

        doid_omim = {}
        for term in onto.get_termobject_list():
            omims = term.get_xrefs('OMIM')
            if omims:
                doid_omim[term.go_id] = omims

        mim_gene = {}
        mim2gene_list = requests.get(MIM2GENE).text.splitlines()

        for line in mim2gene_list:  # Loop from Dima @ Princeton
            if line.startswith('#'):
                continue
            toks = line.split('\t')
            if len(toks) < 3:
                logger.error('Can\'t parse line: %s', line)
                continue
            mim, gtype, gid = toks[:3]
            if gtype in LIMIT_TYPE:
                if mim in mim_gene:
                    logger.warning("MIM already exists: %s", mim)
                if gid:
                    mim_gene[mim] = gid

        mimdiseases = {}
        genemap_list = requests.get(GENEMAP).text.splitlines()
        genemap_version = None

        # TODO: Add support for publications
        for l in genemap_list:  # Loop from Dima @ Princeton
            if l.startswith('#'):
                if l.startswith('# Generated:'):
                    genemap_version = l.split('# Generated:')[1].strip()
                continue

            l_split = l.split('\t')
            status = l_split[6].strip()
            mim_geneid = l_split[8].strip()
            disorders = l_split[11].strip()

            if disorders != '' and status in LIMIT_STATUS and mim_geneid in mim_gene:
                logger.info('%s with status %s', disorders, status)

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
                            if mim_phetype == LIMIT_PHENO:
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

        entrez_gid = {}
        # Iterate over doid->omim mapping
        for doid, omim_list in doid_omim.iteritems():
            term = onto.get_term(doid)
            if term is None:
                continue

            logger.debug("Processing %s", term)

            for omim_id in omim_list:
                if omim_id not in mimdiseases:
                    continue
                mim_entry = mimdiseases[omim_id]
                d_or_s = 'S' if mim_entry.is_susceptibility else 'D'

                for g in mim_entry.genetuples:
                    entrez = g[0]
                    term.add_annotation(gid=entrez, ref=None)

        self._onto = onto
        return onto

if __name__ == '__main__':
    omim = OMIM()
    onto = omim.load()

    onto.print_to_gmt_file('test.txt')
