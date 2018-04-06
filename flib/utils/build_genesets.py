import argparse
import logging

from flib.core.entrez import Entrez
from flib.core.hgmd import HGMD
from flib.core.omim import OMIM
from flib.core.gwas import GWASCatalog
from flib.core.goa import GOA
from flib.core.onto import DiseaseOntology

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(
    description='Generate a file of updated disease gene annotations')
parser.add_argument('--output', '-o', dest='output', type=str,
                                help='output file')
parser.add_argument('--propagate', '-p', dest='propagate', action='store_true',
                    default=False,
                    help='propagate annotations')
parser.add_argument('--namespace', '-n', dest='namespace',
                    default=None,
                    help='term namespace')
parser.add_argument('--databases', '-d', dest='databases',
                    choices=['HGMD', 'OMIM', 'GWASCAT', 'GO'],
                    default=['GO'],
                    nargs='*',
                    help='list of databases')

parser.add_argument('--evidence', '-e', dest='evidence',
                    choices=['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'],
                    default=None,
                    nargs='*',
                    help='list of evidence codes')


args = parser.parse_args()

dbs = set(args.databases)

# Load entrez gene mapping
entrez_map = Entrez()
entrez_map.load()


if 'GO' in dbs:
    goa = GOA()
    onto = goa.load_onto(idmap=entrez_map.get_xref_map())
else:
    onto = DiseaseOntology.generate()
    if 'HGMD' in dbs:
        # Load HGMD annotations
        hgmd = HGMD(host='127.0.0.1', port=3308, user='root', passwd='hgmd')
        hgmd.load_onto(idmap=entrez_map.get_symbol_map(), onto=onto)
    if 'OMIM' in dbs:
        # Load OMIM annotations
        omim = OMIM()
        omim.load_onto(onto=onto)
    if 'GWASCAT' in dbs:
        gwas = GWASCatalog()
        gwas.load_onto(idmap=entrez_map.get_symbol_map(), onto=onto)

if args.evidence:
    logger.info('Filtering by evidence: %s', args.evidence)
    onto.filter_annotations(args.evidence)

if args.propagate:
    logger.info('Propagating gene annotations')
    onto.propagate()

onto.print_to_gmt_file(args.output, p_namespace=args.namespace)
