import argparse
from flib.core.entrez import Entrez
from flib.core.hgmd import HGMD
from flib.core.omim import OMIM
from flib.core.gwas import GWASCatalog
from flib.core.onto import DiseaseOntology

parser = argparse.ArgumentParser(description='Generate a file of updated disease gene annotations')
parser.add_argument('--output', '-o', dest='output', type=str,
                                help='output file')
parser.add_argument('--propagate', '-p', dest='propagate', action='store_true',
                                default=False,
                                help='propagate annotations')
parser.add_argument('--databases', '-d', dest='databases', action='store_true',
                                choices=['HGMD','OMIM','GWASCAT'],
                                default=['HGMD','OMIM'],
                                nargs='*',
                                help='list of disease databases')

args = parser.parse_args()

dbs = set(args.databases)

# Load entrez gene mapping
entrez_map = Entrez()
entrez_map.load()

do = DiseaseOntology.generate()

if 'HGMD' in dbs:
    # Load HGMD annotations
    hgmd = HGMD(host='127.0.0.1', port=3308, user='root', passwd='hgmd')
    hgmd.load_onto(idmap=entrez_map.get_symbol_map(), onto = do)
if 'OMIM' in dbs:
    # Load OMIM annotations
    omim = OMIM()
    omim.load_onto(onto=do)
if 'GWASCAT' in dbs:
    gwas = GWASCatalog()
    gwas.load_onto(idmap=entrez_map.get_symbol_map(), onto = do)

if args.propagate:
    do.propagate()

do.print_to_gmt_file(args.output)
