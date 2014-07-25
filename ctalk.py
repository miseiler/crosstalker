import sys, auc, auc_perm, scripts, clustio, treeio, getopt, random, os

import numpy   as N
import cPickle as cp

from pyvisml import VisML
from random import sample
from itertools import combinations as comb


ITER_ENH  = 100  # Number of permutations to run for enhancement calculation. Note that higher numbers will vastly increase processing time.
ITER_PERM = 1000 # Number of 4/500 permutations to run for significance calculation


def calculate_sa(s, fn):

    thresh = len(s) * 0.85
    nd = dict([ (s.gene_names[i], int(s.M.sum(0)[i] >= thresh)) for i in xrange(len(s.gene_names)) ])
    clustio.write_table(nd, 'gene_presence/%s_top85_gt_1.txt' % fn)

def get_sa(fn):

    sa = clustio.read_table('gene_presence/%s_top85_gt_1.txt' % fn)                
    for x in sa:                    
        sa[x] = int(sa[x])                          
    return sa

def calculate_auc(fn, fln, pathway_dict, path_names):

    sa = get_sa(fn)
    c = clustio.parsers.NullParser()                                             
    c.gene_names = N.array(path_names)       
    c.samples = [ clustio.parsers.SampleData(sample_id=x) for x in c.gene_names ] 
    Q = auc.mp_auc_matrix(fln, [ pathway_dict[x] for x in path_names ], sa, similarity=True)           
    c.M = Q        
    clustio.write_normal(c, 'auc_results/%s_results_reweight_RAW.txt' % fn)

def calculate_perm_test(fn, fln):

    # TODO Hide output
    
    sa = get_sa(fn)
    results = []                                                                 
    for i in xrange(ITER_PERM):
        p1 = sample(fln.gene_names, 4)
        p2 = sample(fln.gene_names, 500)                         
        Q  = auc.mp_auc_matrix(fln, [p1, p2], sa, similarity=True, procs=1)
        results.append(N.sqrt(Q[0][1] * Q[1][0]))                
    clustio.write_list(results, 'auc_results/%s_perm_test_4_500.txt' % fn)

def calculate_sig_connections(fn):

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

def calculate_enhancement(fn1, fn2, pathway_dict, path_lengths):

    # XXX Calculated threshold is more appropriate here to reduce time spent
    q1, q2 = auc_perm.permcomp(fn1, fn2, fln, pathway_dict, path_lengths, threshold=0.0, iter=ITER_ENH)
    M = N.abs(q1 - q2)
    c = clustio.parsers.NullParser()
    c.gene_names = N.array([ str(x) for x in lens ])
    c.samples = [ clustio.parsers.SampleData(sample_id=x) for x in c.gene_names ]
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
    path_names   = list(set(pathway_dict.keys()))

    print('Generating pathway length list...')
    path_lengths = list(set([ len(pathway_dict[x]) for x in pathway_dict ])) 

    fn1 = os.path.splitext(os.path.basename(settings['test']))[0]
    fn2 = os.path.splitext(os.path.basename(settings['control']))[0]

    # Final things needed are reweight_RAW, sig_connections, test_vs_control_thresh_95.txt
    # Steps:
    # calculate_sa creates gene_presence/%s_top85_gt_1.txt
    # * calculate_auc creates auc_results/%s_results_reweight_RAW.txt
    # * calculate_perm_test creates auc_results/%s_perm_test_4_500.txt
    # @ calculate_sig_connections creates sig_connections/%s_sig_connections_95.txt
    # calculate_enhancement creates auc_results/%s_vs_%s_thresh_95.txt
    
    # * denotes activities that can occur simultaneously
    # @ denotes I have tested this to make sure it works

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
                print('New pathway definitions detected! Attempting to remove old results...')
                for f in ('auc_results/%s_results_reweight_RAW.txt' % fn, 'auc_results/%s_perm_test_4_500.txt' % fn, 'sig_connections/%s_sig_connections_95.txt' % fn, 'auc_results/%s_vs_%s_thresh_95.txt' % (fn1, fn2)):
                    try:
                        os.remove(f)
                    except:
                        pass

        if os.path.exists('auc_results/%s_results_reweight_RAW.txt' % fn):
            print('AUC results not found, building...')
            calculate_auc(fn, fln, pathway_dict, path_names)
        else:
            print('Found AUC results')

        if not os.path.exists('auc_results/%s_perm_test_4_500.txt' % fn):
            print('Permutation test results not found, building...')
            calculate_perm_test(fn, fln)
        else:
            print('Found permutation test results')
    
        if not os.path.exists('sig_connections/%s_sig_connections_95.txt' % fn):
            print('Significant connection matrix not found, building...')
            calculate_sig_connections(fn)
        else:
            print('Found significant connection matrix')

    if not os.path.exists('auc_results/%s_vs_%s_thresh_95.txt' % (fn1, fn2)):
        print('Enhancement permutation test results not found, building...')
        calculate_enhancement(fn1, fn2, pathway_dict, path_lengths)
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

    # TODO: AUC and perm test can be done simultaneously

    settings = handle_opts()
    dirstruct()
    fn1, fn2, pathway_dict = generate_missing(settings)

    treeio.create_tree(fn1, fn2, pathway_dict)

    print('done!')
