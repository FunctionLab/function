from obo import OBO

DO_URL = 'https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid-non-classified.obo'
DO_NAME = 'Disease Ontology'

GO_URL = 'http://geneontology.org/ontology/go.obo'
GO_NAME = 'Gene Ontology'


class Ontology:

    @staticmethod
    def generate(obo_file=None, obo_url=None):
        onto = OBO()
        if obo_file:
            onto.load_obo(file=obo_file)
        elif obo_url:
            onto.load_obo(obo_url, remote_location=True, timeout=5)
        return onto


class DiseaseOntology():

    @staticmethod
    def generate():
        return Ontology.generate(obo_url=DO_URL)


class GeneOntology:

    @staticmethod
    def generate():
        return Ontology.generate(obo_url=GO_URL)
