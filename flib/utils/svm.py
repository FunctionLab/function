import argparse

import logging

from multiprocessing import Pool

from flib.core.dab import Dab
from flib.core.gmt import GMT
from flib.core.onto import DiseaseOntology, GeneOntology
from flib.core.labels import OntoLabels, Labels
from flib.core.svm import NetworkSVM

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(
    description='Generate a file of updated disease gene annotations')
parser.add_argument('--input', '-i', dest='input', type=str,
                    required=True,
                    help='Input dab file')
parser.add_argument('--output', '-o', dest='output', type=str,
                                required=True,
                                help='Output directory')
parser.add_argument('--gmt', '-g', dest='gmt', type=str,
                    help='Input GMT (geneset) file')
parser.add_argument('--dir', '-d', dest='dir', type=str,
                    help='Directory of labels')
parser.add_argument('--all', '-a', dest='predict_all', action='store_false',
                    default=True,
                    help='Predict all genes')
parser.add_argument('--slim', '-s', dest='slim', type=str,
                    help='File of slim terms')
parser.add_argument('--threads', '-t', dest='threads', type=int,
                    default=12,
                    help='Number of threads')
parser.add_argument('--best-params', '-b', dest='best_params',
                    action='store_true',
                    default=False,
                    help='Select best parameters by cross validation')
parser.add_argument('--ontology', '-y', dest='ontology',
                    choices=['GO', 'DO'],
                    default='DO',
                    help='Ontology to use for propagation')
parser.add_argument('--namespace', '-n', dest='namespace',
                    help='Limit predictions to terms in the \
                    specified namespace')
parser.add_argument('--flat-output', '-f', dest='flat',
                    action='store_true',
                    default=False,
                    help='Flatten output')

args = parser.parse_args()

MIN_POS, MAX_POS = 5, 300

onto = None
if args.ontology == 'DO':
    logger.info('Loading Disease Ontology')
    onto = DiseaseOntology.generate()
elif args.ontology == 'GO':
    logger.info('Loading Gene Ontology')
    onto = GeneOntology.generate()

if args.gmt:
    # Load GMT genes onto Disease Ontology and propagate
    gmt = GMT(filename=args.gmt)

    # Filter terms by geneset size
    terms = [termid for termid, genes in gmt.genesets.items()
             if len(genes) >= MIN_POS and len(genes) <= MAX_POS]

    # Filter terms by namespace
    if args.namespace and onto:
        for termid in terms[:]:
            term = onto.get_term(termid)
            if term and term.namespace != args.namespace:
                terms.remove(termid)
                logger.info('Ignoring term %s with namespace %s',
                            termid, term.namespace)

    logger.info('Total terms: %i', len(terms))

    if args.slim and onto:
        # Build ontology aware labels
        onto.populate_annotations_from_gmt(gmt)

        lines = open(args.slim).readlines()
        slim_terms = set([l.strip() for l in lines])
        labels = OntoLabels(obo=onto, slim_terms=slim_terms)
    else:
        labels = Labels(gmt=gmt)

elif args.dir:
    labels = Labels(labels_dir=args.dir)
    terms = [term for term in labels.get_terms()
             if len(labels.get_labels(term)[0]) >= MIN_POS and
             len(labels.get_labels(term)[0]) <= MAX_POS]
else:
    logger.error('Insufficient options to proceed. \
            Please provide a GMT file or a directory of labels')
    exit()

dab = Dab(args.input)
svm = NetworkSVM(dab, preload=args.predict_all)


def run_svm(term):
    (pos, neg) = labels.get_labels(term)

    logger.info('Running SVM for %s, %i pos, %i neg', term, len(pos), len(neg))

    svm.predict(pos, neg,
                predict_all=args.predict_all,
                best_params=args.best_params)
    svm.print_predictions(args.output + '/' + term,
                          pos, neg, term, flat=args.flat)


if args.threads > 1:
    pool = Pool(args.threads)
    pool.map(run_svm, terms)
    pool.close()
else:
    for term in terms:
        run_svm(term)
