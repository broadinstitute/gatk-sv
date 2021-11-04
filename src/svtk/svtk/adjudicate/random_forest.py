#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve


def rf_classify(metrics, trainable, testable, features, labeler, cutoffs, name,
                clean_cutoffs=False):
    """Wrapper to run random forest and assign probabilities"""
    rf = RandomForest(trainable, testable, features, cutoffs, labeler, name,
                      clean_cutoffs)
    rf.run()
    metrics.loc[rf.testable.index, name] = rf.probs
    cutoffs = rf.cutoffs.copy()

    # evidence = name.split('_')[0]

    #  rf.clean.to_csv('{0}_training.txt'.format(evidence), index=False, sep='\t')
    del rf.clean
    del rf.testable
    del rf.rf
    del rf

    return cutoffs


class RandomForest:
    def __init__(self, trainable, testable, features, cutoffs, labeler, name,
                 clean_cutoffs=False, max_train_size=100000):
        def has_null_features(df):
            return df[features].isnull().any(axis=1)

        self.clean = trainable.loc[~has_null_features(trainable)].copy()
        if self.clean.shape[0] == 0:
            raise Exception('No clean variants found')

        self.testable = testable.loc[~has_null_features(testable)].copy()

        self.features = features

        self.labeler = labeler
        self.encoder = LabelEncoder().fit(['Fail', 'Pass'])

        self.name = name
        self.clean_cutoffs = clean_cutoffs
        self.cutoff_features = cutoffs
        self.cutoffs = None
        self.max_train_size = max_train_size

    def run(self):
        sys.stderr.write('Labeling training data...\n')
        self.label_training_data()
        sys.stderr.write('Selecting training data...\n')
        self.select_training_data()
        sys.stderr.write('Learning probabilities...\n')
        self.learn_probs()
        sys.stderr.write('Learning cutoffs...\n')
        self.learn_cutoffs()
        sys.stderr.write('Trimming probabilities...\n')
        self.cutoff_probs()

    def label_training_data(self):
        self.clean['label'] = self.labeler.label(self.clean)

    def select_training_data(self):
        self.train = self.clean.loc[self.clean.label != 'Unlabeled']

        if self.train.shape[0] >= self.max_train_size:
            max_subset_size = int(self.max_train_size / 2)

            passes = self.train.loc[self.train.label == 'Pass']
            if passes.shape[0] >= max_subset_size:
                passes = passes.sample(max_subset_size)

            fails = self.train.loc[self.train.label == 'Fail']
            if fails.shape[0] >= max_subset_size:
                fails = fails.sample(max_subset_size)

            self.train = pd.concat([passes, fails])

        if self.train.loc[self.train.label == 'Pass'].shape[0] == 0:
            raise Exception('No Pass variants included in training set')
        if self.train.loc[self.train.label == 'Fail'].shape[0] == 0:
            raise Exception('No Fail variants included in training set')

    def learn_probs(self):
        X_train = self.train[self.features].values

        y_train = self.encoder.transform(self.train.label)

        self.rf = RandomForestClassifier(n_estimators=500, random_state=343124,
                                         oob_score=True, max_features=None)

        self.rf.fit(X_train, y_train)

        X = self.testable[self.features].values
        probs = self.rf.predict_proba(X)

        self.probs = probs[:, 1]

    def learn_cutoffs(self):
        cutoffs = {}
        # Restrict learning cutoffs to "clean" variants
        if self.clean_cutoffs:
            cutoff_metrics = self.clean
        else:
            cutoff_metrics = self.testable

        cutoff_metrics['pass_cutoffs'] = True
        passing = cutoff_metrics['pass_cutoffs']

        for feature in self.cutoff_features['indep']:
            metric = cutoff_metrics[feature]
            idx = np.searchsorted(self.testable.index, cutoff_metrics.index)
            if self.name == "PE_prob":
                cutoff1 = learn_cutoff_dist(metric, self.probs[idx])
                cutoff2 = learn_cutoff_fdr(metric, self.probs[idx])
                cutoff = min([cutoff1, cutoff2])
            else:
                cutoff = learn_cutoff_dist(metric, self.probs[idx])

            cutoffs[feature] = cutoff
            passing = passing & (cutoff_metrics[feature] >= cutoff)

        passing = cutoff_metrics.loc[passing]
        for feature in self.cutoff_features['dep']:
            metric = passing[feature]
            # Subset probabilities to those in passing set
            idx = np.searchsorted(self.testable.index, passing.index)
            if self.name == "PE_prob":
                cutoff1 = learn_cutoff_dist(metric, self.probs[idx])
                cutoff2 = learn_cutoff_fdr(metric, self.probs[idx])
                cutoffs[feature] = min([cutoff1, cutoff2])
            else:
                cutoffs[feature] = learn_cutoff_dist(metric, self.probs[idx])

        self.cutoffs = pd.DataFrame.from_dict({'cutoff': cutoffs},
                                              orient='columns')\
            .reset_index()
        self.cutoffs = self.cutoffs.rename(columns=dict(index='metric'))

    def cutoff_probs(self):
        self.testable['prob'] = self.probs
        passes = self.testable.prob >= 0.5

        self.testable['passes_all_cutoffs'] = True

        # If metrics are below the observed cutoff, force failure
        for idx, row in self.cutoffs.iterrows():
            metric, cutoff = row['metric'], row['cutoff']
            self.testable.loc[passes & (
                self.testable[metric] < cutoff), 'prob'] = 0.499

            self.testable['passes_all_cutoffs'] = self.testable.passes_all_cutoffs & (
                self.testable[metric] >= cutoff)

        self.testable.loc[self.testable.passes_all_cutoffs &
                          (self.testable.prob < 0.5), 'prob'] = 0.501

        self.probs = self.testable.prob.values


def learn_cutoff_dist(metric, probs):
    preds = metric.values

    # Pass/fail if greater/less than 0.5
    classify = np.vectorize(lambda x: 1 if x >= 0.5 else 0)
    truth = classify(probs)

    # If all variants which passed prior cutoffs also passed random forest,
    # return minimum value instead of trying to compute cutoff
    if 0 not in truth:
        return preds.min()

    fpr, tpr, thresh = roc_curve(truth, preds)
    dist = np.sqrt((fpr - 0) ** 2 + (tpr - 1) ** 2)
    best_idx = np.argmin(dist)

    # If cutoff set at no instances, scikit-learn sets thresh[0] to max(y_score) + 1
    if best_idx == 0:
        return thresh[best_idx] - 1

    return thresh[best_idx]


def learn_cutoff_fdr(metric, probs, fdr_cff=.05):
    preds = metric.values

    # Pass/fail if greater/less than 0.5
    classify = np.vectorize(lambda x: 1 if x >= 0.5 else 0)
    truth = classify(probs)

    # If all variants which passed prior cutoffs also passed random forest,
    # return minimum value instead of trying to compute cutoff
    if 0 not in truth:
        return preds.min()

    fpr, tpr, thresh = roc_curve(truth, preds)

    for i in range(len(fpr)):
        if fpr[i] > fdr_cff:
            break
    best_idx = i
    # If cutoff set at no instances, scikit-learn sets thresh[0] to max(y_score) + 1
    if best_idx == 0:
        return thresh[best_idx] - 1

    return thresh[best_idx]
