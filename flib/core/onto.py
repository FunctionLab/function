from obo import OBO

from flib import settings

class Ontology:

    @staticmethod
    def generate(obo_file=None, obo_url=None):
        onto = OBO()
        if obo_file:
            onto.load_obo(file=obo_file)
        elif obo_url:
            onto.load_obo(obo_url, remote_location=True, timeout=15)
        return onto


class DiseaseOntology():

    @staticmethod
    def generate():
        return Ontology.generate(obo_url=settings.DO_URL)


class GeneOntology:

    @staticmethod
    def generate():
        return Ontology.generate(obo_url=settings.GO_URL)
