import argparse
import numpy as np
import random
import os
from collections import namedtuple
from operator import itemgetter

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict, GridSearchCV, KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score, classification_report, average_precision_score, brier_score_loss
from sklearn.calibration import CalibratedClassifierCV, _SigmoidCalibration
from sklearn.isotonic import IsotonicRegression

from dab import Dab
from gmt import GMT
from omim import OMIM
from onto import DiseaseOntology
from labels import OntoLabels

class NetworkSVM:

    default_params = {'C':50}
    tuned_parameters = [
        {'C': [.0001, .001, .01, .1, 1, 10, 100], 'class_weight':['balanced', None]},
    ]

    def __init__(self, dab):
        self._dab = dab

    def _dab_matrix(self):
        if not self._X_all:
            # Load dab as matrix
            self._X_all = np.empty([self._dab.get_size(), self._dab.get_size()])
            for i, g in enumerate(self._dab.gene_list):
                if not i % 1000:
                    logger.info('Loaded %i', i)
                self._X_all[i] = self._dab.get(g)
        return self._X_all

    def predict(self, pos_genes, neg_genes,
            predict_all=False,
            best_params=False,
            prob_fit='SIGMOID',
            cv_folds=5):

        logger.info("Running %i fold SVM", cv_folds)

        # Group training genes
        train_genes = [g for g in (pos_genes | neg_genes) if self._dab.get_index(g) is not None]
        train_genes_idx = [self._dab.get_index(g) for g in train_genes]

        # Subset training matrix and labels
        if predict_all:
            X = self._dab_matrix()[train_genes_idx]
            y = np.array([1 if g in pos_genes else -1 for g in train_genes])
        else:
            X = np.empty([len(train_genes), self._dab.get_size()])
            y = np.empty(len(train_genes))
            for i, g in enumerate(train_genes):
                X[i] = self._dab.get(g)
                y[i] = 1 if g in pos_genes else -1

        params = NetworkSVM.default_params

        if best_params:
            # Set the parameters by cross-validation
            score = 'average_precision'
            clf = GridSearchCV(LinearSVC(), NetworkSVM.tuned_parameters, cv=3, n_jobs=10,
                               scoring=score)
            clf.fit(X, y)
            params = clf.best_params_

        train_scores, train_probs = np.empty(len(train_genes)), np.empty(len(train_genes))
        train_scores[:], train_probs[:] = np.NAN, np.NAN
        scores, probs = None, None

        kf = StratifiedKFold(n_splits=cv_folds)
        for cv, (train, test) in enumerate(kf.split(X, y)):
            X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]

            logger.info('Learning SVM')
            clf = LinearSVC(**params)
            clf.fit(X_train, y_train)

            logger.info('Predicting SVM')
            if predict_all:
                scores_cv = clf.decision_function(X_all)
                scores = scores_cv if scores is None else np.column_stack((scores, scores_cv))

                for idx in test:
                    train_scores[idx] = scores_cv[train_genes_idx[idx]]
            else:
                scores_cv = clf.decision_function(X_test)
                for i,idx in enumerate(test):
                    train_scores[idx] = scores_cv[i]

        if prob_fit == 'ISO':
            ir = IsotonicRegression(out_of_bounds='clip')
            Y = label_binarize(y, [-1,1])
            ir.fit(train_scores, Y[:,0])
            train_probs = ir.predict(train_scores)
        else:
            Y = label_binarize(y, [-1,1])
            sc = _SigmoidCalibration()
            sc.fit(train_scores, Y)
            train_probs = sc.predict(train_scores)

        if predict_all:
            scores = np.median(scores, axis=1)
            for i,idx in enumerate(train_genes_idx):
                scores[idx] = train_scores[i]

            probs = np.median(probs, axis=1)
            for i,idx in enumerate(train_genes_idx):
                probs[idx] = train_probs[i]

            genes = dab.gene_list
        else:
            scores = train_scores
            genes = train_genes
            probs = train_probs


        self._predictions = sorted(zip(genes, scores, probs), key=itemgetter(1), reverse=True)

        return self._predictions

    def print_predictions(self, ofile, pos_genes, neg_genes):
        with open(ofile, 'w') as outfile:
            for (g,s,p) in self._predictions:
                if g in pos_genes:
                    label = '1'
                elif g in neg_genes:
                    label = '-1'
                else:
                    label = '0'
                line = [g, label, str(s), str(p), '\n']
                outfile.write('\t'.join(line))
            outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a file of updated disease gene annotations')
    parser.add_argument('--input', '-i', dest='input', type=str,
                                    help='input dab file')
    parser.add_argument('--output', '-o', dest='output', type=str,
                                    help='output directory')
    parser.add_argument('--gmt', '-g', dest='gmt', type=str,
                                    help='input GMT (geneset) file')
    parser.add_argument('--dir', '-d', dest='dir', type=str,
                                    help='directory of labels')
    parser.add_argument('--all', '-a', dest='all', action='store_true',
                                    default=False,
                                    help='predict all genes')
    parser.add_argument('--best-params', '-b', dest='best_params', action='store_true',
                                    default=False,
                                    help='select best parameters by cross validation')
    parser.add_argument('--prob', '-p', dest='prob_fit',
                                    choices=['SIGMOID','ISO'],
                                    default=None,
                                    help='probability fit')
    args = parser.parse_args()

    dab = Dab(args.input)
    svm = NetworkSVM(dab)

    do = DiseaseOntology.generate()

    if args.gmt:
        gmt = GMT(filename=args.gmt)
        for (gsid, genes) in gmt.genesets.iteritems():
            term = do.get_term(gsid)
            for gid in genes:
                term.add_annotation(gsid, gid)
    else:
        OMIM().load_onto(onto=do)

    do.propagate()

    lines = open('../../files/do_slim.txt').readlines()
    slim_terms = set([l.strip() for l in lines])

    ol = OntoLabels(obo=do, slim_terms=slim_terms)

    for term in do.get_termobject_list():
        if len(term.annotations) == 0:
            continue

        (pos, neg) = ol.get_labels(term.go_id)
        if len(pos) < 5:
            continue
        print term.go_id, len(pos), len(neg)
        predictions = svm.predict(pos, neg)

