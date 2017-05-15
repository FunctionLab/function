import sys
import os
import math
import numpy as np
import string
from optparse import OptionParser
from collections import defaultdict

from sklearn.metrics import roc_auc_score, average_precision_score

usage = "usage: %prog [options]"
parser = OptionParser(usage, version = "%prog dev-unreleased")
parser.add_option("-d", "--dir", dest="dir", help="svmperf directory", metavar="FILE")

(options, args) = parser.parse_args()

files = os.listdir(options.dir)
files.sort()

for f in files:
    labels, scores = [], []

    for l in open(options.dir + '/' + f):
        gene, status, val, prob = l.strip().split('\t')[0:4]
        val = prob
        if status != '0':
            if status == '1':
                labels.append(True)
            elif status == '-1':
                labels.append(False)
            scores.append(float(val))

    labels, scores = np.array(labels), np.array(scores)

    print f, average_precision_score(labels, scores), roc_auc_score(labels, scores)
