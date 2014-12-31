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
import sys
import clustio, scripts

import numpy as N
import multiprocessing as mp
from itertools import combinations as comb
from itertools import combinations_with_replacement as combr

from auc import roc, mutual, auc
from ctalk import get_sa, mutualize

MP_MAX_QUEUE = 16


def findpathwaysizes(fn1, fn2, pathway_dict, sizes, threshold):
    mat1 = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % fn1)
    mat2 = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % fn2)
    assert (mat1.gene_names == mat2.gene_names).all()
    mutualize(mat1)
    mutualize(mat2)
    M = N.abs(mat1.M - mat2.M)
    res = []
    nd = {}
    for i in xrange(len(sizes)):
        nd[sizes[i]] = i
    for i, j in comb(xrange(len(M)), 2):
        if M[i][j] >= threshold:
            res.append(list(sorted([ nd[len(pathway_dict[mat1.gene_names[i]])], nd[len(pathway_dict[mat1.gene_names[j]])] ])))
    #return list(set([ tuple(x) for x in res if (N.array(x) > 16).all() and (N.array(x) < 350).all() ]))
    return list(set([ tuple(x) for x in res ]))

def permcomp_by_size(fn1, fn2, fln, p1size, p2size, iter=1000):
    sa1 = get_sa(fn1)
    sa2 = get_sa(fn2)
    assert set(sa1.keys()) == set(sa2.keys())
    results = []
    seedmat = N.random.rand(2, 2)
    Q1 = mp_auc_matrix(fln, [p1size, p2size], sa1, similarity=True, seedmat = seedmat, iter=iter)
    Q2 = mp_auc_matrix(fln, [p1size, p2size], sa2, similarity=True, seedmat = seedmat, iter=iter)
    return Q1, Q2

def permcomp(fn1, fn2, fln, pathway_dict, sizes, threshold=0.0, procs=mp.cpu_count(), iter=100):
    sa1 = get_sa(fn1)
    sa2 = get_sa(fn2)
    assert set(sa1.keys()) == set(sa2.keys())

    # Dec 30 2014 currently on hold while we work out how this is going to go
    #pairs = findpathwaysizes(fn1, fn2, pathway_dict, sizes, threshold)
    #print('Calculating permutations for %s/%s possible pairs' % (len(pairs), (len(sizes) * (len(sizes) + 1)) / 2))
    
    seedmat = N.random.rand(len(sizes))
    result1 = mp_auc_matrix(fln, sizes, sa1, similarity=True, seedmat=seedmat, procs=procs, iter=iter)
    result2 = mp_auc_matrix(fln, sizes, sa2, similarity=True, seedmat=seedmat, procs=procs, iter=iter)
    return result1, result2

def predictability_perm_roc(s, size1, size2, overlap, sa, iter, similarity, seed):
    """
    
    Calculate ROC of predictability for pathway sizes size1 and size2
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

    import random
    random.seed(seed)
    from random import sample

    res = []
    for _ in xrange(iter):

        # Take a random sample without replacement of total necessary size
        gene_pool = sample(c.gene_names, size1 + size2 - overlap)

        # Find splits
        div1 = size1 - overlap
        div2 = size2 - overlap + div1

        # Engineer two sets of size size1 and size2 with exactly overlap genes overlapping
        gl1 = gene_pool[:div1]     + gene_pool[div2:]
        gl2 = gene_pool[div1:div2] + gene_pool[div2:]
    
        #Sens = TP / (TP + FN)
        #Spec = FP / (TN + FP)
        #Plot sens vs 1 - spec

        roc1 = roc(c, gl1, gl2, sa)
        roc2 = roc(c, gl2, gl1, sa)

        res.append(mutual(roc1, roc2))

    return res

def _pr(q, rq, ns, num, sa):

    #print('Worker started')

    res = []

    while True:
        
        v = q.get()
        if v is None:
            #print('Worker received termination request! Stopping...')
            break

        psize1, psize2, overlap = v
        #if i % 100 == 0 and j == i+1:
        #    print('Current job status: %s' % i)
    
        seed = ns.seed_dict[(psize1, psize2, overlap)]

        perms = predictability_perm_roc(ns.s, psize1, psize2, overlap, sa, ns.iter, ns.similarity, seed)
        
        res.append((psize1, psize2, overlap, perms))
        #res.append((psize2, psize1, overlap, perms))

    #print('Worker %s writing to file...' % num)
    #clustio.write_list(['\t'.join([str(elem) for elem in x]) for x in res], 'results_%s.txt' % num)
    
    #print('Queuing result (%s)' % num)
    rq.put(res)

def mp_auc_matrix(s, pathwaysizes, sa, similarity=False, iter=1000, procs=mp.cpu_count(), seed_dict=None):
    """
    Expects a list of lists of pathways (groups of elements) found in s.gene_names
    Returns an asymmetric matrix of AUC values where M[i][j] is the predictive value of pathway i for pathway j
    
    """
    
    if seed_dict is None:
        print('Random seeds for paired perm tests not given, generating new random seeds')
        seed_dict = dict(zip(pathwaysizes, N.random.rand(len(pathwaysizes))))
    
    assert len(seed_dict) == len(pathwaysizes)

    print('Performing permutation tests for %s sets of pathway sizes and overlaps' % len(pathwaysizes))

    result_dict = {}

    ns = mp.Manager()
    q  = mp.Queue()
    rq = mp.Queue()
    
    ns.s = s
    ns.similarity = similarity
    ns.iter = iter
    ns.seed_dict = seed_dict
    
    # Set up queue

    qsplit = [ pathwaysizes[MP_MAX_QUEUE*i:MP_MAX_QUEUE*(i+1)] for i in xrange(len(pathwaysizes)/MP_MAX_QUEUE + 1) ] # Splits the queue up into sizes MP_MAX_QUEUE, plus a remainder list. I don't know why it works.

    for combs in qsplit:

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
            
            for i, j, k, v in res:
                result_dict[(i,j,k)] = v
                result_dict[(j,i,k)] = v

    return result_dict
