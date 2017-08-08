import sys
import os
import math
import numpy as np
import string
from optparse import OptionParser
from collections import defaultdict

from sklearn.metrics import roc_auc_score, average_precision_score

usage = "usage: %prog [options]"
parser = OptionParser(usage, version="%prog dev-unreleased")
parser.add_option("-d", "--dir", dest="dir",
    help="svmperf directory", metavar="FILE")
parser.add_option("-g", "--neg-genes", dest="neg_genes",
    help="gene list of negative examples", metavar="FILE")
parser.add_option("-l", "--label-col", dest="label_col",
    default=1,
    type=int,
    help="column of the label (zero-indexed)")
parser.add_option("-s", "--score-col", dest="score_col",
    default=2,
    type=int,
    help="column of the score (zero-indexed)")



(options, args) = parser.parse_args()

neg_genes = None
if options.neg_genes:
    neg_genes = set()
    with open(options.neg_genes) as f:
        for g in f.readlines():
            neg_genes.add(g.strip())

files = os.listdir(options.dir)
files.sort()

for f in files:
    labels, scores, probs = [], [], []

    for l in open(options.dir + '/' + f):
        tok = l.strip().split('\t')
        gene, label, score = tok[0], tok[options.label_col], tok[options.score_col]
        if label != '0' or neg_genes is not None:
            if label == '1':
                labels.append(True)
                scores.append(float(score))
            elif label == '-1' or (neg_genes and gene in neg_genes):
                labels.append(False)
                scores.append(float(score))

    labels, scores, probs = np.array(labels), np.array(scores), np.array(probs)

    print f, len([l for l in labels if l]), len([l for l in labels if not l]), \
        average_precision_score(labels, scores), roc_auc_score(labels, scores)
