
from go import go

DO_URL = 'https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid-non-classified.obo'
DO_NAME = 'Disease Ontology'


class DiseaseOntology:

    @staticmethod
    def generate(obo_file=None):
        do = go()
        if obo_file:
            do.load_obo(file=obo_file)
        else:
            do.load_obo(DO_URL, remote_location=True, timeout=5)
        return do
