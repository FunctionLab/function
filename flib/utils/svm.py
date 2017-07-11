import argparse

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from multiprocessing import Pool

from flib.core.dab import Dab
from flib.core.gmt import GMT
from flib.core.omim import OMIM
from flib.core.onto import Ontology, DiseaseOntology, GeneOntology
from flib.core.labels import OntoLabels, Labels
from flib.core.svm import NetworkSVM

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
parser.add_argument('--all', '-a', dest='predict_all', action='store_true',
                    default=False,
                    help='Predict all genes')
parser.add_argument('--slim', '-s', dest='slim', type=str,
                    help='File of slim terms')
parser.add_argument('--threads', '-t', dest='threads', type=int,
                    default=12,
                    help='Number of threads')
parser.add_argument('--best-params', '-b', dest='best_params', action='store_true',
                    default=False,
                    help='Select best parameters by cross validation')
parser.add_argument('--ontology', '-y', dest='ontology',
                    choices=['GO', 'DO'],
                    default='DO',
                    nargs=1,
                    help='Ontology to use for propagation')
args = parser.parse_args()

MIN_POS, MAX_POS = 5, 500

if args.ontology == 'DO':
    onto = DiseaseOntology.generate()
elif args.ontology == 'GO':
    onto = GeneOntology.generate()
else:
    onto = Ontology.generate()

if args.gmt:
    # Load GMT genes onto Disease Ontology and propagate
    gmt = GMT(filename=args.gmt)
    onto.populate_annotations_from_gmt(gmt)
    onto.propagate()

    # Filter terms by number of gene annotations
    terms = [term.go_id for term in onto.get_termobject_list()
             if len(term.annotations) >= MIN_POS and len(term.annotations) <= MAX_POS]

    if args.slim:
        # Build ontology aware labels
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
svm = NetworkSVM(dab)


def run_svm(term):
    (pos, neg) = labels.get_labels(term)

    logger.info('Running SVM for %s, %i pos, %i neg', term, len(pos), len(neg))

    predictions = svm.predict(pos, neg,
                              predict_all=args.predict_all,
                              best_params=args.best_params)
    svm.print_predictions(args.output + '/' + term, pos, neg)

pool = Pool(args.threads)
pool.map(run_svm, terms)
pool.close()
