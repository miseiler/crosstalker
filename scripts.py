"""

Copyright 2009 Michael Seiler
Rutgers University
miseiler@gmail.com

This file is part of ConsensusCluster.

ConsensusCluster is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ConsensusCluster is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ConsensusCluster.  If not, see <http://www.gnu.org/licenses/>.


"""

import os
import numpy as N

from itertools import combinations as comb
from clustio.parsers import SampleData
from copy import deepcopy

from clustio import *


def transpose_sdata(s):

    s.M = s.M.T.astype(N.float32)

    gl   = list(s.gene_names)
    sids = list([ x.sample_id for x in s.samples ])

    s.samples    = [ parsers.SampleData(sample_id = x) for x in gl ]
    s.gene_names = N.array(sids)

    return s

def scale_to_set(sdata, *filenames):
    """

    scale_to_set(filename)

        Removes all but those sample_ids you specifiy.

        filenames    - filenames or dicts
                       each file containing list of sample ids to use
                       or each dict containing name->list of sample ids

    Returns: modified sdata object, dict of cluster->indices

    """

    new_sdata = parsers.NullParser()

    defined_clusters = list_or_files(*filenames)

    sample_id_list    = [ x.sample_id for x in sdata.samples ]
    samples_to_keep   = sum([ defined_clusters[x] for x in defined_clusters ], [])
    sample_indices    = argintersect(sample_id_list, samples_to_keep)
    
    sample_classes    = dict([ (defined_clusters[k][i], k) for k in defined_clusters
                               for i in xrange(len(defined_clusters[k])) ])

    #Adjustment
    new_sdata.samples = [ SampleData(cluster_id=sdata.samples[i].cluster_id, sample_id=sdata.samples[i].sample_id,
                                     sample_num=sdata.samples[i].sample_num, index=sdata.samples[i].index,
                                     sample_class=sample_classes[sdata.samples[i].sample_id]) for i in sample_indices ]

    new_sdata.M = sdata.M.take(tuple(sample_indices), 0)
    new_sdata.gene_names = sdata.gene_names.copy()

    sample_id_list = new_sdata.sample_ids
    
    for name in defined_clusters: #If samples aren't in the main, ignore them
        sample_list = defined_clusters[name]
        def_indices = argintersect(sample_list, sample_id_list)
        defined_clusters[name] = [ sample_list[x] for x in def_indices ]

    return new_sdata, defined_clusters

def scale_probes(sdata, *filenames):
    """

    scale_probes(sdata, filename)
        
        Removes all gene probes except those you specify

        filename    - File(s) containing a list of probes, one on each line
                      Also accepted: Lists, dicts.  Only the values will be used in the dicts.

    Returns: modified sdata object

    """

    new_sdata = parsers.NullParser()

    plist = sum(list_or_files(*filenames).values(), [])

    probes_to_keep = tuple(argintersect(sdata.gene_names, plist))

    new_sdata.M = sdata.M.take(probes_to_keep, 1)
    new_sdata.gene_names = sdata.gene_names.take(probes_to_keep)
    new_sdata.samples = [ deepcopy(x) for x in sdata.samples ]

    return new_sdata

def new_defined_clusters(sdata, conv):
    """

    Define different clusters than the ones specified by your Defined Clusters, whether
    through the GUI, modification of keep_list, or through the command line

    sdata: sample data obj
    conv: conversion dict, keys sample ids values new cluster assignments
    Stick this in your preprocess function (see common.py for subclassing help)

    """
    
    new_clusts = {}
    s_ids = [x.sample_id for x in sdata.samples]

    for s_id in s_ids:
        if s_id in conv:
            new_clusts.setdefault(conv[s_id], []).append(s_id)
        else:
            new_clusts.setdefault('Unknown', []).append(s_id)

    return new_clusts

def create_network_from_sdata(s, weighted=False):

    import networkx as nx

    # Are we dealing with a similarity graph?
    assert s.M.shape[0] == s.M.shape[1]
    assert (N.array([ x.sample_id for x in s.samples ]) == s.gene_names).all()

    G = nx.Graph()

    G.add_nodes_from(s.gene_names)

    for i, j in comb(xrange(len(s)), 2):
        weight = s.M[i][j]

        if weight:
            if weighted:
                G.add_edge(s.gene_names[i], s.gene_names[j], weight=weight)
            else:
                G.add_edge(s.gene_names[i], s.gene_names[j])

    return G

def create_visml_from_sdata(s, distance='euclidean', dmatrix=False, allow_zero_weight=False):

    from pyvisml import VisML
    
    if not dmatrix:
        import scipy
        dm = eval('scipy.spatial.distance.%s' % distance)
    
        # Distance matrix M
        M = N.zeros((s.M.shape[0], s.M.shape[0]), dtype=N.float32)
        for i, j in comb(xrange(len(s)), 2):
            M[i][j] = dm(s.M[i], s.M[j])
    
        # Normed to [0,1]
        M /= M.max()
    else:
        M = s.M

    tree = VisML.create_empty_visml_tree(layout='elegant:100', fineArt='False')

    nodes = []
    for name in s.sample_ids:
        nodes.append(tree.add_node(name, '0', '0'))

    for i, j in comb(xrange(len(s)), 2):
        if M[i][j] or allow_zero_weight:
            tree.add_edge(nodes[i], nodes[j], weight=str(M[i][j]))

    return tree

def create_sdata_from_visml(tree, methods=[]):
    """
    Creates a symmetric SampleData object from tree, incorporating only those links in methods if given

    """

    from pyvisml import VisML

    A = VisML.VisMLTree(tree)
    
    Anodes = {}
    for k in A.nodes:
        Anodes[k.name] = k

    sample_ids = Anodes.keys()
    sample_ids.sort()

    M = N.zeros((len(sample_ids), len(sample_ids)), dtype=N.float32)

    Midx = {}
    for i in xrange(len(sample_ids)):
        Midx[sample_ids[i]] = i

    for k in sample_ids:
        for link in Anodes[k].links:
            
            if methods:
                if link.method not in methods:
                    continue

            i = Midx[k]
            j = Midx[link.to]

            if link.weight is not None:
                M[i][j] = float(link.weight)
            else:
                M[i][j] = 1.0

    assert (M == M.T).all()

    c = parsers.NullParser()
    c.M = M
    c.gene_names = N.array(sample_ids)
    c.samples = [ parsers.SampleData(sample_id=x) for x in c.gene_names ]

    return c
