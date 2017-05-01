import sys 
import os
import MySQLdb as mdb 
from do import DiseaseOntology
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

EVIDENCE = frozenset(['DM','DFP'])

class HGMD:

    def __init__(self, host='127.0.0.1', port=3306, user=None, passwd=None):
        self._host = host
        self._port = port
        self._user = user
        self._passwd = passwd
        self._data = None

    def load_data(self): 
        data = set()
        try:
            con = mdb.connect(host=self._host, port=self._port, 
                    user=self._user, passwd=self._passwd)
            cur = con.cursor()

            cur.execute("select gene, tag, phenotype, cui from hgmd_pro.allmut join hgmd_phenbase.hgmd_mutation on hgmd_pro.allmut.acc_num = hgmd_phenbase.hgmd_mutation.acc_num join hgmd_phenbase.hgmd_phenotype on hgmd_phenbase.hgmd_phenotype.phen_id = hgmd_phenbase.hgmd_mutation.phen_id join hgmd_phenbase.phenotype_concept on hgmd_phenbase.hgmd_phenotype.phen_id = hgmd_phenbase.phenotype_concept.phen_id")
            rows = cur.fetchall()
            for row in rows:
                #(acc_num,gene_sym,phen_id1,phen_id2,rela,cui,phen_id3,phenotype) = row
                gene,tag,phenotype,cui = row 
                result = (gene, cui, phenotype,tag)
                data.add(result)
            con.close()

        except mdb.Error, e:
            logger.error('Error quering HGMD database')
            return False 

        self._data = data

        return True 

    def load(self, onto=None, evidence=EVIDENCE):
        if not self._data:
            self.load_data()

        if not onto:
            onto = DiseaseOntology.generate()

        xrefs = onto.get_xref_mapping('UMLS_CUI')

        for (gene, cui, phenotype, evd) in self._data:
            if cui not in xrefs or evd not in evidence:
                continue

            for doid in xrefs[cui]:
                term = onto.get_term(doid)
                term.add_annotation(gid=gene)

        self._onto = onto

        return onto

if __name__ == '__main__':
    hgmd = HGMD(host='127.0.0.1', port=3308, user='root', passwd='hgmd')
    onto = hgmd.load()
    onto.print_to_gmt_file('test.txt')
