"""
Entrez Gene
"""
GENEINFO_URLS = {
    'Homo sapiens':
    'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
}

UNIPROT_URLS = {
    'Homo sapiens':
    'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/' + \
        'knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'
}

UNIPROT_PREFIX = 'UniProtKB'


"""
GOA - Gene Ontology Annotations
"""
GOA_NAMES = {
    'Arabidopsis thaliana': 'tair',
    'Homo sapiens': 'human',
    'Mus musculus': 'mgi',
    'Rattus norvegicus': 'rgd',
    'Danio rerio': 'zfin',
    'Drosophila melanogaster': 'fb',
    'Saccharomyces cerevisiae': 'sgd',
    'Caenorhabditis elegans': 'wb',
    'Pseudomonas aeruginosa': 'pseudocap'
}

GOA_PREFIX = ['goa_', 'gene_association.']
GOA_ASSOC_SUFFIX = ['.gaf.gz', '.gz']
GOA_INFO_SUFFIX = ['.gaf.json', '.json']

GOA_ASSOC_URL = 'http://www.geneontology.org/gene-associations/'
GOA_VERSION_KEY = 'submissionDate'


"""
GWAS Catalog
"""
GWAS_URL = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'


"""
HGMD
"""
HGMD_DEFAULT_EVIDENCE = frozenset(['DM', 'DFP'])


"""
OMIM - Online Mendelian Inheritance in Man
"""
OMIM_KEY = 'DUBJEOOeQLWrZuLTgSmSyA'
OMIM_MIM2GENE = 'http://omim.org/static/omim/data/mim2gene.txt'
OMIM_GENEMAP = 'http://data.omim.org/downloads/' + OMIM_KEY + '/genemap.txt'

OMIM_LIMIT_TYPE = set(['gene', 'gene/phenotype'])
OMIM_LIMIT_PHENO = '(3)'
OMIM_LIMIT_STATUS = ['C', 'P']

"""
DO - Disease Ontology
"""
DO_URL = 'https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid-non-classified.obo'
DO_NAME = 'Disease Ontology'


"""
GO - Gene Ontology
"""
GO_URL = 'http://geneontology.org/ontology/go.obo'
GO_NAME = 'Gene Ontology'
