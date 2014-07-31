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
import sys, auc, auc_perm, scripts, clustio, treeio, getopt, random, os

import numpy   as N
import cPickle as cp

from pyvisml import VisML
from random import sample
from itertools import combinations as comb


ITER_ENH  = 100  # Number of permutations to run for enhancement calculation. Note that higher numbers will vastly increase processing time.
ITER_PERM = 100  # Number of permutations to run for significance calculation
ITER_ENH_T= 1000 # Number of permutations to run to find a minimum threshold for calculating enhancement for that size pair


def _calculate_perm_test(fn, fln):
    # OLD METHOD, SUPERCEDED ON JULY 29 2014

    # TODO Hide output
    
    sa = get_sa(fn)
    results = []                                                                 
    for i in xrange(ITER_PERM):
        p1 = sample(fln.gene_names, 4)
        p2 = sample(fln.gene_names, 500)                         
        Q  = auc.mp_auc_matrix(fln, [p1, p2], sa, similarity=True, procs=1)
        results.append(N.sqrt(Q[0][1] * Q[1][0]))                
    clustio.write_list(results, 'auc_results/%s_perm_test_4_500.txt' % fn)

def _calculate_sig_connections(fn):
    # OLD METHOD, SUPERCEDED ON JULY 29 2014

    s = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % fn)
    for i, j in comb(xrange(len(s)), 2):
        v = N.sqrt(s.M[i][j] * s.M[j][i])
        s.M[i][j] = s.M[j][i] = v
    results = [ float(x) for x in clustio.read_list('auc_results/%s_perm_test_4_500.txt' % fn) ]
    results.sort()
    t95 = results[int(0.95 * ITER_PERM)]
    for i, j in comb(xrange(len(s)), 2):
        if s.M[i][j] < t95:
            s.M[i][j] = s.M[j][i] = 0
    clustio.write_normal(s, 'sig_connections/%s_sig_connections_95.txt' % fn)


def get_sa(fn):

    sa = clustio.read_table('gene_presence/%s_top85_gt_1.txt' % fn)                
    for x in sa:                    
        sa[x] = int(sa[x])                          
    return sa

def create_new_symm_dset(features):

    c = clustio.parsers.NullParser()                                             
    c.gene_names = N.array(features)       
    c.samples = [ clustio.parsers.SampleData(sample_id=x) for x in c.gene_names ] 

    return c

def mutualize(s):

    for i, j in comb(xrange(len(s)), 2):
        v = N.sqrt(s.M[i][j] * s.M[j][i])
        s.M[i][j] = s.M[j][i] = v

def calculate_sa(s, fn):

    M = (s.M >= 1.0).sum(0)
    thresh = len(s) * 0.85
    nd = dict([ (s.gene_names[i], int(M[i] >= thresh)) for i in xrange(len(s.gene_names)) ])
    clustio.write_table(nd, 'gene_presence/%s_top85_gt_1.txt' % fn)

def calculate_auc(fn, fln, pathway_dict, path_names):

    sa = get_sa(fn)
    c = create_new_symm_dset(path_names)
    Q = auc.mp_auc_matrix(fln, [ pathway_dict[x] for x in path_names ], sa, similarity=True)           
    c.M = Q        
    clustio.write_normal(c, 'auc_results/%s_results_reweight_RAW.txt' % fn)

def calculate_perm_test(fn, fln, path_lengths):

    sa = get_sa(fn)
    M  = auc_perm.mp_auc_matrix(fln, path_lengths, sa, similarity=True, iter=ITER_PERM)
    M.sort(2)
    c = create_new_symm_dset(path_lengths)
    c.M = M[:,:,int(0.95 * ITER_PERM)]
    clustio.write_normal(c, 'auc_results/%s_perm_test.txt' % fn)

def calculate_sig_connections(fn, pathway_dict):

    s = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % fn)
    mutualize(s)
    c = clustio.ParseNormal('auc_results/%s_perm_test.txt' % fn)
    
    psizes = dict([ (x, len(pathway_dict[x])) for x in pathway_dict ])
    cidx   = dict([ (int(i), c.gene_names[i]) for i in xrange(len(c.gene_names)) ])

    for i, j in comb(xrange(len(s)), 2):
        ilen = psizes[s.gene_names[i]]
        jlen = psizes[s.gene_names[j]]
        thresh = c.M[cidx[ilen]][cidx[jlen]]
        if s.M[i][j] < thresh:
            s.M[i][j] = s.M[j][i] = 0

    clustio.write_normal(s, 'sig_connections/%s_sig_connections_95.txt' % fn)

def calculate_enhancement_threshold(fn1, fn2, fln, path_lengths):

    pl = sorted(path_lengths, reverse=True)
    q1, q2 = auc_perm.permcomp_by_size(fn1, fn2, fln, pl[0], pl[1], iter=ITER_ENH_T)
    r2 = N.abs(q1[0][1] - q2[0][1])
    r2.sort()
    return r2[int(0.95 * ITER_ENH_T)]

def calculate_enhancement(fn1, fn2, fln, pathway_dict, path_lengths, threshold=0.0):

    q1, q2 = auc_perm.permcomp(fn1, fn2, fln, pathway_dict, path_lengths, threshold=threshold, iter=ITER_ENH)
    M = N.abs(q1 - q2)
    c = create_new_symm_dset(path_lengths)
    M.sort(2)
    c.M = M[:,:,int(0.95 * ITER_ENH)]
    clustio.write_normal(c, 'auc_results/%s_vs_%s_thresh_95.txt' % (fn1, fn2))
    clustio.write_normal(c, 'auc_results/%s_vs_%s_thresh_95.txt' % (fn2, fn1))

def dirstruct():
    """
    Build directory structure:
    auc_results
    gene_presence
    sig_connections
    
    """

    for dr in ('auc_results', 'gene_presence', 'sig_connections'):
        if not os.path.isdir(dr):
            os.mkdir(dr)

def handle_opts():

    # TODO: Defaults go here!
    settings = {'test': None, 'control': None, 'FLN': None, 'pathways': None}

    def usage(err=None):
        print('\nUSAGE: python ctalk.py [OPTIONS]\n')
        print('\t-f, --fln <filename>\t\t\tUse the FLN at <filename>, rather than the default. Tab-delimited with symmetric row and col headers.')
        print('\t-p, --pathways <filename>\t\tUse pathway definitions at <filename>, rather than the default (Oct 2013 KEGG definitions). Should be a pickled python dictionary.')
        print('\t-t, --test_condition <filename>\t\tTest condition data file. Tab-delimited with samples on cols and features on rows, with headers on both. e.g., tumor data.')
        print('\t-c, --control_condition <filename>\tControl condition data file. Same format as --test_condition. e.g., normal tissue data.')
        print('\n\tEXAMPLE: python ctalk.py --pathways hsa_paths --fln FLN_hsa.txt --test_condition LumA_tcga_data.txt --control_condition Normal_tcga_data.txt\n')

        if err is not None:
            print(err)
            print

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hf:p:t:c:', ['fln=', 'pathways=', 'test_condition=', 'control_condition=', 'help'])
    except getopt.GetoptError as err:
        usage(err)
        sys.exit(2)

    for o, a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif o in ('-f', '--fln'):
            settings['FLN'] = a
        elif o in ('-p', '--pathways'):
            settings['pathways'] = a
        elif o in ('-t', '--test_condition'):
            settings['test'] = a
        elif o in ('-c', '--control_condition'):
            settings['control'] = a
        else:
            usage('Option not recognized: %s' % o)

    if None in settings.values():
        usage()
        print('Missing data: %s' % ' '.join([ x for x in settings if settings[x] is None ]))
        sys.exit(1)

    return settings

def generate_missing(settings):
    """
    Figure out what data files are missing towards our final network and fill them in

    """

    print('Loading FLN...')
    fln = clustio.ParseNormal(settings['FLN'])

    print('Loading pathway definitions...')
    f = open(settings['pathways'], 'r')
    pathway_dict = cp.load(f)
    f.close()
    
    print('Found %s pathways' % len(pathway_dict))
    path_names   = sorted(pathway_dict.keys())

    print('Generating pathway length list...')
    path_lengths = sorted(set([ len(pathway_dict[x]) for x in pathway_dict ])) 

    fn1 = os.path.splitext(os.path.basename(settings['test']))[0]
    fn2 = os.path.splitext(os.path.basename(settings['control']))[0]

    # Final things needed are reweight_RAW, sig_connections, test_vs_control_thresh_95.txt
    # Steps:
    # @ calculate_sa creates gene_presence/%s_top85_gt_1.txt
    # @ * calculate_auc creates auc_results/%s_results_reweight_RAW.txt
    # * calculate_perm_test creates auc_results/%s_perm_test.txt
    # @ calculate_sig_connections creates sig_connections/%s_sig_connections_95.txt
    # @ calculate_enhancement creates auc_results/%s_vs_%s_thresh_95.txt
    
    # * denotes activities that can occur simultaneously
    # @ denotes I have tested this to make sure it works in a manner consistent to earlier code

    for fn, datafile in [(fn1, settings['test']), (fn2, settings['control'])]:
        print
        print('Looking for required %s files and generating them if needed...' % fn)
        print
    
        if not os.path.exists('gene_presence/%s_top85_gt_1.txt' % fn):
            print('Gene presence definitions not found, building...')
            calculate_sa(clustio.ParseNormal(datafile), fn)
        else:
            print('Found gene presence definitions')
    
        if os.path.exists('auc_results/%s_results_reweight_RAW.txt' % fn):
            s = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % fn)
            if not s.sample_ids == path_names:
                del s
                print('New pathway definitions detected! Attempting to remove old results...')
                for f in ('auc_results/%s_results_reweight_RAW.txt' % fn, 'auc_results/%s_perm_test_4_500.txt' % fn, 'sig_connections/%s_sig_connections_95.txt' % fn, 'auc_results/%s_vs_%s_thresh_95.txt' % (fn1, fn2)):
                    try:
                        os.remove(f)
                    except:
                        pass

        if not os.path.exists('auc_results/%s_results_reweight_RAW.txt' % fn):
            print('AUC results not found, building...')
            calculate_auc(fn, fln, pathway_dict, path_names)
        else:
            print('Found AUC results')

        if not os.path.exists('auc_results/%s_perm_test.txt' % fn):
            print('Permutation test results not found, building...')
            calculate_perm_test(fn, fln, path_lengths)
        else:
            print('Found permutation test results')
    
        if not os.path.exists('sig_connections/%s_sig_connections_95.txt' % fn):
            print('Significant connection matrix not found, building...')
            calculate_sig_connections(fn, pathway_dict)
        else:
            print('Found significant connection matrix')

    if not os.path.exists('auc_results/%s_vs_%s_thresh_95.txt' % (fn1, fn2)):
        print('Enhancement permutation test results not found')
        print('Calculating an appropriate threshold to reduce permutation search space...')
        threshold = calculate_enhancement_threshold(fn1, fn2, fln, path_lengths)
        print('Building enhancement permutation test matrix...')
        calculate_enhancement(fn1, fn2, fln, pathway_dict, path_lengths, threshold=threshold)
    else:
        print('Found enhancement permutation test results')

    return fn1, fn2, pathway_dict


if __name__ == '__main__':

    """

    Complete crosstalk workflow:
    
    REQUIRES FROM USER:

    pathway definitions
    tumor conditions
    normal condition
    FLN
    threshold?

    DO:
    create sig_connections, auc_results, gene_presence
    accept txt file, convert to sdata
    do 4 - 500 permutation test
    do AUC calculations
    do tumor - normal AUC calculations
    build networks

    As we work, look for files that already exist and announce them. Don't do that work if found.

    """

    # TODO: By sharing a seed matrix, we can do just one permutation test for both enh and AUC threshold

    settings = handle_opts()
    dirstruct()
    fn1, fn2, pathway_dict = generate_missing(settings)

    treeio.create_tree(fn1, fn2, pathway_dict)

    print('done!')
