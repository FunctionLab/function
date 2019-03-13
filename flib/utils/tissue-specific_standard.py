import argparse

import logging
from functools import reduce
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from collections import defaultdict

from flib.core.obo import OBO
from flib.core.gmt import GMT
from flib.core.dab import Dab

"""
Based on Arjun's perl script:
https://github.com/FunctionLab/function-perl/blob/master/scripts/construct_tspfc-undir_c1234_gold-std.pl
"""
parser = argparse.ArgumentParser(
    description='Constructs gold-standard across tissues given a template of global gold-std \
    edges and tissue-specific gene-expression.')
parser.add_argument('--positives', '-p', dest='pos',
                    #required=True,
                    help='Global gold-std positive edges in DAT format')
parser.add_argument('--negatives', '-n', dest='neg',
                    #required=True,
                    help='Global gold-std negative edges in DAT format')
parser.add_argument('--global', '-g', dest='global_std',
                    #required=True,
                    help='Global gold-std edges in DAB format')
parser.add_argument('--tissue-genes', '-t', dest='tissue_genes',
                    #required=True,
                    help='Direct tissue-specific gene-expression annotations in GMT format')
parser.add_argument('--ubiq-genes', '-u', dest='ubiq_genes',
                    help='List of ubiquitous genes')
parser.add_argument('--tissue-onto', '-o', dest='tissue_onto',
                    help='Tissue ontology (e.g. Brenda)')
parser.add_argument('--background', '-b', dest='back_genes',
                    help='List of background genes')

args = parser.parse_args()

ubiq = None
if args.ubiq_genes:
    ubiq = set()
    with open(args.ubiq_genes) as f:
        for l in f.readlines():
            ubiq.add(l.strip())

    logger.info('Total ubiquitous genes: %i', len(ubiq))

onto = None
if args.tissue_onto:
    onto = OBO(args.tissue_onto)

tissue_genes = GMT(args.tissue_genes)

if onto:
    onto.populate_annotations_from_gmt(tissue_genes)
    onto.propagate()
    tissue_genes = onto.as_gmt()

with open(args.pos) as f:
    edge_lines = f.readlines()

tissue_std_edges = defaultdict(dict)

for line in edge_lines:
    g1, g2, std = line.strip().split()[:3]
    edge = frozenset([g1, g2])

    # Skip edges where both genes are ubiquitous
    if ubiq and len(edge & ubiq) == 2:
        continue

    for (tissue, genes) in tissue_genes.genesets.items():
        # Include positive edges where both genes are expressed in the tissue
        # Include if one gene is tissue-specific and one gene is ubiq
        if len(edge & genes) == 2 or \
            ubiq and len(edge & genes) == 1 and list((edge - genes))[0] in ubiq:

            if std not in tissue_std_edges[tissue]:
                tissue_std_edges[tissue][std] = set()

            tissue_std_edges[tissue][std].add(edge)

within_edges_pos = reduce(set.union, [tissue_std_edges[tissue]['1'] \
    for tissue in tissue_std_edges])
within_edges_neg = reduce(set.union, [tissue_std_edges[tissue]['0'] \
    for tissue in tissue_std_edges])
tissue_all_pos = reduce(frozenset.union, within_edges_pos)

for (tissue, stds) in tissue_std_edges.items():
    c3 = tissue_std_edges[tissue]['0']
    for c3_edge in list(c3):
        # Limit negative edges to positive genes
        if len(c3_edge & tissue_all_pos) < 2:
            c3.remove(c3_edge)

for (tissue, stds) in tissue_std_edges.items():
    c1 = tissue_std_edges[tissue]['1']
    c1_g = reduce(frozenset.union, c1)
    c1_tspc = c1_g - ubiq

    c3 = tissue_std_edges[tissue]['0']
    c3_g = reduce(frozenset.union, c3)
    c3_tspc = c3_g - ubiq

    c2 = within_edges_pos - c1
    c4 = within_edges_neg - c3
    if onto:
        for dtissue in onto.get_descendents(tissue):
            if dtissue in tissue_std_edges:
                c2 -= tissue_std_edges[dtissue]['1']
                c4 -= tissue_std_edges[dtissue]['0']
        for atissue in onto.get_ancestors(tissue):
            if atissue in tissue_std_edges:
                c2 -= tissue_std_edges[atissue]['1']
                c4 -= tissue_std_edges[atissue]['0']

    print(tissue, len(c1_tspc), len(c1), len(c2), len(c3), len(c4))

