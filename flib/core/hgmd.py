import sys
import os
import MySQLdb as mdb
from onto import DiseaseOntology
from entrez import Entrez
import logging

from flib.settings import HGMD_DEFAULT_EVIDENCE

logging.basicConfig()
logger = logging.getLogger(__name__)

class HGMD:

    def __init__(self, host='127.0.0.1', port=3306, user=None, passwd=None):
        self._host = host
        self._port = port
        self._user = user
        self._passwd = passwd
        self._data = None
        self._onto = None

    def load_data(self):
        data = set()
        try:
            con = mdb.connect(host=self._host, port=self._port,
                              user=self._user, passwd=self._passwd)
            cur = con.cursor()

            cur.execute('select gene, tag, phenotype, cui from hgmd_pro.allmut '
                        'join hgmd_phenbase.hgmd_mutation on '
                        'hgmd_pro.allmut.acc_num = hgmd_phenbase.hgmd_mutation.acc_num '
                        'join hgmd_phenbase.hgmd_phenotype on '
                        'hgmd_phenbase.hgmd_phenotype.phen_id = hgmd_phenbase.hgmd_mutation.phen_id '
                        'join hgmd_phenbase.phenotype_concept on '
                        'hgmd_phenbase.hgmd_phenotype.phen_id = hgmd_phenbase.phenotype_concept.phen_id')
            rows = cur.fetchall()
            for row in rows:
                gene, tag, phenotype, cui = row
                result = (gene, cui, phenotype, tag)
                data.add(result)
            con.close()

        except mdb.Error as e:
            logger.error('Error quering HGMD database')
            return False

        self._data = data

        return True

    def load_onto(self, onto=None, evidence=HGMD_DEFAULT_EVIDENCE, idmap=None):
        if not self._data:
            self.load_data()

        if not onto:
            onto = DiseaseOntology.generate()

        xrefs = onto.get_xref_mapping('UMLS_CUI')

        for (gene, cui, phenotype, evd) in self._data:
            if cui not in xrefs or evd not in evidence:
                continue

            mapped_genes = idmap[gene] if idmap else (gene,)

            for doid in xrefs[cui]:
                term = onto.get_term(doid)
                for gid in mapped_genes:
                    term.add_annotation(gid=gid)

        self._onto = onto

        return onto

if __name__ == '__main__':
    entrez_map = Entrez()
    entrez_map.load()

    hgmd = HGMD(host='127.0.0.1', port=3308, user='root', passwd='hgmd')
    onto = hgmd.load_onto(idmap=entrez_map.get_symbol_map())

    onto.print_to_gmt_file('test.txt')
