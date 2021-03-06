"""

Copyright 2014 Michael Seiler
Boston University
miseiler@gmail.com

This file is part of Crosstalker.

Crosstalker is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Crosstalker is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Crosstalker.  If not, see <http://www.gnu.org/licenses/>.


"""

# Jul 30 2014 deprecated CROC method
# import croc
import sys
import clustio, scripts

import numpy as N
import multiprocessing as mp
from itertools import combinations as comb

try:
    from numpy import isclose
except:
    def isclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        """
        Returns a boolean array where two arrays are element-wise equal within
        a tolerance.
    
        This function was added to numpy v1.7.0, and the version you are
        running has been backported from numpy v1.8.1. See its documentation
        for more details.
        """
        def within_tol(x, y, atol, rtol):
            with N.errstate(invalid='ignore'):
                result = N.less_equal(abs(x-y), atol + rtol * abs(y))
            if N.isscalar(a) and N.isscalar(b):
                result = bool(result)
            return result
    
        x = N.array(a, copy=False, subok=True, ndmin=1)
        y = N.array(b, copy=False, subok=True, ndmin=1)
        xfin = N.isfinite(x)
        yfin = N.isfinite(y)
        if all(xfin) and all(yfin):
            return within_tol(x, y, atol, rtol)
        else:
            finite = xfin & yfin
            cond = N.zeros_like(finite, subok=True)
            # Since we're using boolean indexing, x & y must be the same shape.
            # Ideally, we'd just do x, y = broadcast_arrays(x, y). It's in
            # lib.stride_tricks, though, so we can't import it here.
            x = x * N.ones_like(cond)
            y = y * N.ones_like(cond)
            # Avoid subtraction with infinite/nan values...
            cond[finite] = within_tol(x[finite], y[finite], atol, rtol)
            # Check for equality of infinite values...
            cond[~finite] = (x[~finite] == y[~finite])
            if equal_nan:
                # Make NaN == NaN
                cond[N.isnan(x) & N.isnan(y)] = True
            return cond

def roc(s, seed_list, target_list, sa):

    """
    Takes a similarity matrix s, a seed gene list, a target gene list,

    Returns a set of (score, label) tuples, where label is 1 for pos and 0 for neg
    Suitable for CROC ScoreData object

    """
    
    #try:
    #    assert s.get_samples(seed_list).M.any()
    #    assert s.get_samples(target_list).M.any()
    #except:
    #    print(seed_list)
    #    print(target_list)
    #    raise ValueError, 'Unable to find seed and/or target lists in sdata!'

    gl1set = set(seed_list)
    gl2set = set(target_list)

    genes  = set(s.gene_names)

    union = gl1set & gl2set

    # Consider only genes which are not in the seed list, unless they are found in both seed and target list
    q = s.get_samples(seed_list).get_features(list( (genes ^ (genes & gl1set)) | union ))
    weights = q.M.sum(0)

    c = 1000000 # Dummy gene starter value
    wadd = []
    gadd = []

    # Reweight genes which overlap and are present in the sample set to be infinite
    # Reweight genes which appear in the target set to be 0, and replace them with dummy targets
    for i in xrange(len(weights)):
        if q.gene_names[i] in union and sa[q.gene_names[i]]:
            weights[i] = sys.maxint
        if q.gene_names[i] in gl2set and not sa[q.gene_names[i]]:
            
            # Add gene to end, with 0 weight
            wadd.append(0)
            gadd.append(q.gene_names[i])

            # Current weight gene becomes dummy
            q.gene_names[i] = '%s' % c
            c += 1

    weights = list(weights) + wadd
    q.gene_names = list(q.gene_names) + gadd

    return N.array(weights), N.array([ int(x in gl2set) for x in q.gene_names ])

    # Jul 30 2014 return changed for new AUC method
    #weights = [ (weights[i], q.gene_names[i] in gl2set) for i in xrange(len(weights)) ]

    #return weights

def predictability_roc(s, gl1, gl2, sa, similarity=False):
    """
    
    Calculate ROC of predictability for gl1 and gl2
    given s, an sdata object that is assumed to be a distance matrix normalized between 0 and 1

    if similarity is True, the matrix is assumed to be a similarity matrix instead

    """

    try:
        assert (s.gene_names == s.sample_ids).all()
    except:
        raise ValueError, 'sdata object is not a distance/similarity matrix'

    try:
        assert not (s.M < 0).any()
        assert not (s.M > 1).any()
    except:
        raise ValueError, 'Unnormalized matrix; data found which is outside [0,1] bound'


    if not similarity:
        Q = (1 - s.M.copy()) # Convert to similarity matrix
    else:
        Q = s.M.copy()

    # Create a new s object for convenience, certainly not efficient
    c = clustio.parsers.NullParser()
    c.samples = [ clustio.parsers.SampleData(sample_id=x) for x in s.gene_names ]
    c.gene_names = s.gene_names.copy()
    c.M = Q

    #Sens = TP / (TP + FN)
    #Spec = FP / (TN + FP)
    #Plot sens vs 1 - spec

    roc1 = roc(c, gl1, gl2, sa)
    roc2 = roc(c, gl2, gl1, sa)

    return roc1, roc2

def mutual(roc1, roc2):
    """

    Calculates the mutual predictability of gl1 and gl2

    """

    return N.sqrt(auc(roc1) * auc(roc2))

def _auc(roc, sweep_method='smooth'):
    """

    Calculate the area under the curve using the CROC library
    Superceded Jul 30 2014

    """
    
    return croc.ROC(croc.ScoredData(roc).sweep_threshold(sweep_method)).area()

def auc(roc):
    """

    AUC score calculation adapted from the implementation in scikits-learn 0.16

    """
    
    y_score, y_true = roc

    # sort scores and corresponding truth values
    desc_score_indices = N.argsort(y_score, kind="mergesort")[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]

    # y_score typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve.
    # We need to use isclose to avoid spurious repeated thresholds
    # stemming from floating point roundoff errors.
    distinct_value_indices = N.where(N.logical_not(isclose(
        N.diff(y_score), 0)))[0]
    threshold_idxs = N.r_[distinct_value_indices, y_true.size - 1]

    # accumulate the true positives with decreasing threshold
    tps = y_true.cumsum()[threshold_idxs]
    fps = 1 + threshold_idxs - tps
    
    fpr = fps / float(fps[-1])
    tpr = tps / float(tps[-1])
    
    return N.trapz(tpr, fpr)

def _pr(q, rq, ns, num, sa):

    #print('Worker started')

    res = []

    while True:
        
        v = q.get()
        if v is None:
            #print('Worker received termination request! Stopping...')
            break

        i, j = v
        #if i % 100 == 0 and j == i+1:
        #    print('Current job status: %s' % i)
    
        p1 = ns.pathways[i]
        p2 = ns.pathways[j]
    
        roc1, roc2 = predictability_roc(ns.s, p1, p2, sa, ns.similarity)
    
        res.append((i, j, auc(roc1)))
        res.append((j, i, auc(roc2)))

    #print('Worker %s writing to file...' % num)
    #clustio.write_list(['\t'.join([str(elem) for elem in x]) for x in res], 'results_%s.txt' % num)
    
    #print('Queuing result (%s)' % num)
    rq.put(res)

def mp_auc_matrix(s, pathways, sa, similarity=False, procs=mp.cpu_count()):
    """
    Expects a list of lists of pathways (groups of elements) found in s.gene_names
    Returns an asymmetric matrix of AUC values where M[i][j] is the predictive value of pathway i for pathway j
    
    """
    
    M = N.zeros((len(pathways), len(pathways)), N.float32)

    ns = mp.Manager()
    q  = mp.Queue()
    rq = mp.Queue()
    
    ns.s = s
    ns.pathways = pathways
    ns.similarity = similarity

    # Set up queue
    combs = list(comb(xrange(len(pathways)), 2))
    for c in combs:
        q.put(c)

    # Start workers

    #print('Worker count: %s' % procs)
    workers = {}
    for i in xrange(procs):
        workers[i] = mp.Process(target=_pr, args=(q, rq, ns, i, sa))
        workers[i].start()
        q.put(None)

    # If we ask the workers to join, the results queue fills and the threads will block while waiting for the queue to unfill

    # Wait for workers
    for k in workers:
        res = rq.get() # Blocks until results queue completes, which is also when the workers terminate
        
        for i, j, v in res:
            M[i][j] = v

    return M
