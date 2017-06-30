import argparse
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

from optparse import OptionParser

from flib.core.obo import OBO
from flib.core.gmt import GMT

parser = argparse.ArgumentParser(description='Generate propagated gene annotation lists from ontology and association files')
parser.add_argument(
    "-o",
    "--obo-file",
    dest="obo",
    required=True,
    help="A obo file",
    metavar="FILE")
parser.add_argument(
    "-a",
    "--association-file",
    dest="ass",
    help="A gene association file",
    metavar="FILE")
parser.add_argument(
    "-b",
    dest="term_col",
    type=int,
    help="The column of the annotations file containing the term identifiers",
    default=4)
parser.add_argument(
    "-g",
    dest="gcol",
    type=int,
    help="The column of the annotations file containing the desired identifiers",
    default=1)
parser.add_argument(
    "-G",
    "--gmt-file",
    dest="gmt",
    help="GMT file of gene associations")
parser.add_argument(
    "-d",
    "--output-prefix",
    dest="opref",
    help="The prefix for output files",
    metavar="string")
parser.add_argument(
    "-f",
    "--output-filename",
    dest="ofile",
    help="If given outputs all go term/gene annotation pairs to this file, file is created in the output prefix directory.",
    metavar="string")
parser.add_argument(
    "-i",
    "--id-file",
    dest="idfile",
    help="File to map existing gene ids to the desired identifiers in the format <gene id>\\t<desired id>\\n",
    metavar="FILE")
parser.add_argument(
    "-p",
    action="store_true",
    dest="propagate",
    default=False,
    help="Propagate gene annotations")
parser.add_argument(
    "-t",
    "--terms-file",
    dest="terms",
    help="File of terms to limit output",
    metavar="FILE")
parser.add_argument(
    "-n",
    "--namespace",
    dest="nspace",
    default="biological_process",
    help="limit the GO term output to the input namespace: (biological_process, cellular_component, molecular_function)",
    metavar="STRING")
parser.add_argument(
    "-A",
    dest="assoc_format",
    action="store_true",
    default=False,
    help="If we are printing to a file (-f), pass this to get a full association file back.")
parser.add_argument(
    "-u",
    dest="pub_filter",
    action="store_true",
    default=False,
    help="Filter annotations from high-throughput publications (>50 annotations)")

args = parser.parse_args()

if args.obo is None:
    sys.stderr.write("--obo file is required.\n")
    sys.exit()
if args.pub_filter and args.nspace is None:
    sys.stderr.write(
        "--When filtering by publication, must provide GO namespace.\n")
    sys.exit()

id_name = None
if args.idfile is not None:
    id_name = IDMap(args.idfile)

gene_ontology = OBO(args.obo)

logger.info('Populating gene associations')
if args.ass:
    gene_ontology.populate_annotations(
        args.ass,
        gene_col=args.gcol,
        term_col=args.term_col)
elif args.gmt:
    gmt = GMT(args.gmt)
    gene_ontology.populate_annotations_from_gmt(gmt)
else:
    sys.stderr.write("--Provide gene annotations from an association file or a GMT file")
    exit()

if args.pub_filter:
    pub_counts = defaultdict(set)
    for (term_id, term) in gene_ontology.go_terms.iteritems():
        if term.namespace != args.nspace:
            continue
        for a in term.annotations:
            pub_counts[a.ref].add((term, a))
    for (ref, annots) in pub_counts.iteritems():
        if len(annots) > 50:
            logger.info(
                'Removing %i annotations from: %s',
                ref,
                len(annots))
            for (term, a) in annots:
                term.remove_annotation(a)

if args.idfile is not None:
    gene_ontology.map_genes(id_name)

if args.propagate:
    logger.info('Propagating gene associations')
    gene_ontology.propagate()

gterms = None
if args.terms:
    f = open(args.terms, 'r')
    gterms = []
    for line in f:
        termid = line.rstrip('\n')
        gterms.append(termid)
    f.close()

if args.ofile:
    gene_ontology.print_to_single_file(
        args.opref +
        '/' +
        args.ofile,
        gterms,
        args.nspace,
        args.assoc_format)
else:
    gene_ontology.print_to_dir(args.opref,
        gterms,
        args.nspace)
