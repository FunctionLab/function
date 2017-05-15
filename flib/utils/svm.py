import argparse
import numpy as np
import random
from collections import namedtuple
from operator import itemgetter

from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict, GridSearchCV, KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, classification_report, average_precision_score, brier_score_loss
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import BaggingClassifier
from sklearn.calibration import CalibratedClassifierCV

from flib.core.dab import Dab
from flib.core.omim import OMIM
from flib.core.hgmd import HGMD
from flib.core.onto import DiseaseOntology
from flib.core.entrez import Entrez
from flib.core.gmt import GMT

parser = argparse.ArgumentParser(description='Generate a file of updated disease gene annotations')
parser.add_argument('--input', '-i', dest='input', type=str,
                                help='input dab file')
parser.add_argument('--output', '-o', dest='output', type=str,
                                help='output directory')
parser.add_argument('--gmt', '-g', dest='gmt', type=str,
                                help='input GMT (geneset) file')
parser.add_argument('--all', '-a', dest='all', action='store_true',
                                default=False,
                                help='predict all genes')
parser.add_argument('--best-params', '-b', dest='best_params', action='store_true',
                                default=False,
                                help='select best parameters by cross validation')
parser.add_argument('--geneset_id', '-G', dest='geneset_id', type=str,
                                help='geneset id')

parser.add_argument('--prob', '-p', dest='prob_fit',
                                choices=['SIGMOID','ISO'],
                                default=None,
                                help='probability fit')


args = parser.parse_args()



standards = {}
Std = namedtuple('Std', ['pos', 'neg'])

gmt = GMT(filename=args.gmt)
if args.geneset_id:
    pos_genes = gmt.get_genes(args.geneset_id)
    neg_genes = gmt.genes - pos_genes
    standards[args.geneset_id] = Std(pos=pos_genes, neg=neg_genes)
else:
    for (gsid, genes) in gmt.genesets.iteritems():
        pos_genes = gmt.get_genes(gsid)
        neg_genes = gmt.genes - pos_genes
        if len(pos_genes) >= 10 and len(pos_genes) <= 1000:
            standards[gsid] = Std(pos=pos_genes, neg=neg_genes)

print len(standards.keys())

dab = Dab(args.input)
if args.all:
    # Load dab as matrix
    X_all = np.empty([dab.get_size(), dab.get_size()])
    for i, g in enumerate(dab.gene_list):
        if not i % 1000:
            print i
        X_all[i] = dab.get(g)


for gsid, std in standards.iteritems():
    print 'Predicting', gsid, len(std.pos), len(std.neg)
    pos_genes, neg_gens = std.pos, std.neg

    # Group training genes
    train_genes = [g for g in (pos_genes | neg_genes) if dab.get_index(g) is not None]
    train_genes_idx = [dab.get_index(g) for g in train_genes]

    if args.all:
        # Subset training matrix and labels
        X = X_all[train_genes_idx]
        y = np.array([1 if g in pos_genes else -1 for g in train_genes])
    else:
        X = np.empty([len(train_genes), dab.get_size()])
        y = np.empty(len(train_genes))
        for i, g in enumerate(train_genes):
            X[i] = dab.get(g)
            y[i] = 1 if g in pos_genes else -1

    if args.best_params:
        # Set the parameters by cross-validation
        tuned_parameters = [
            {'C': [.0001, .001, .01, .1, 1, 10, 100], 'class_weight':['balanced', None]},
        ]
        score = 'average_precision'

        print("# Tuning hyper-parameters for %s" % score)

        clf = GridSearchCV(LinearSVC(), tuned_parameters, cv=3, n_jobs=10,
                           scoring=score)
        clf.fit(X, y)
        best_params = clf.best_params_

        print(clf.best_params_)
        means = clf.cv_results_['mean_test_score']
        stds = clf.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, clf.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))
    else:
        best_params = {'C':50, 'class_weight':'balanced'}


    train_scores, train_probs = np.empty(len(train_genes)), np.empty(len(train_genes))
    train_scores[:], train_probs[:] = np.NAN, np.NAN
    scores, probs = None, None

    kf = StratifiedKFold(n_splits=5)
    for cv, (train, test) in enumerate(kf.split(X, y)):
        X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]

        print "Learning SVM"
        clf = LinearSVC(**best_params)
        clf.fit(X_train, y_train)

        if args.prob_fit == 'SIGMOID':
            clf_prob = CalibratedClassifierCV(clf, cv=3, method='sigmoid')
            clf_prob.fit(X_train, y_train)
        elif args.prob_fit == 'ISO':
            clf_prob = CalibratedClassifierCV(clf, cv=3, method='isotonic')
            clf_prob.fit(X_train, y_train)

        print "Predicting SVM"
        if args.all:
            scores_cv = clf.decision_function(X_all)
            scores = scores_cv if scores is None else np.column_stack((scores, scores_cv))

            for idx in test:
                train_scores[idx] = scores_cv[train_genes_idx[idx]]

            # Store fitted probabilities
            if args.prob_fit:
                probs_cv = clf_prob.predict_proba(X_all)[:,1]
                probs = probs_cv if probs is None else np.column_stack((probs, probs_cv))

                for idx in test:
                    train_probs[idx] = probs_cv[train_genes_idx[idx]]

        else:
            scores_cv = clf.decision_function(X_test)
            for i,idx in enumerate(test):
                train_scores[idx] = scores_cv[i]

            print roc_auc_score(y_test, scores_cv)

            # Store fitted probabilities
            if args.prob_fit:
                probs_cv = clf_prob.predict_proba(X_test)[:,1]

                print roc_auc_score(y_test, probs_cv)

                for i,idx in enumerate(test):
                    train_probs[idx] = probs_cv[i]
                sorted_scores = sorted(zip([train_genes[i] for i in test], scores_cv, probs_cv), key=itemgetter(1), reverse=True)
                '''
                with open(args.output + '/' + gsid + '_' + str(cv), 'w') as outfile:
                    for g,s,p in sorted_scores:
                        if g in pos_genes:
                            label = '1'
                        elif g in neg_genes:
                            label = '-1'
                        else:
                            label = '0'
                        line = [g, label, str(s), str(p), '\n']
                        outfile.write('\t'.join(line))
                    outfile.close()
                '''

    if args.all:
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

    print 'Performance:', \
        len(neg_genes & set(train_genes)), \
        roc_auc_score(y, train_scores), \
        roc_auc_score(y, train_probs), \
        average_precision_score(y, train_scores), \
        average_precision_score(y, train_probs)

    if args.output:
        sorted_scores = sorted(zip(genes, scores, probs), key=itemgetter(1), reverse=True)
        with open(args.output + '/' + gsid, 'w') as outfile:
            for (g,s,p) in sorted_scores:
                if g in pos_genes:
                    label = '1'
                elif g in neg_genes:
                    label = '-1'
                else:
                    label = '0'
                line = [g, label, str(s), str(p), '\n']
                outfile.write('\t'.join(line))
            outfile.close()
