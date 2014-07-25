from pyvisml import VisML
import clustio, scripts
from itertools import combinations as comb

import numpy as N

NO_COMPARISON_SUFFIX = 'no_comparison'
COMPARISON_SUFFIX    = 'final'

# XXX requires pathway definitions

def difftree(pathway_dict, test_cond_name, control_cond_name, infile_suffix=NO_COMPARISON_SUFFIX):            
    # To be performed on freshly created scripts.create_visml_from_sdata(%s_sig_connections.txt, dmatrix=True)
    
    A = VisML.VisMLTree('%s_%s.xml' % (test_cond_name, infile_suffix))
    B = clustio.ParseNormal('sig_connections/%s_sig_connections_95.txt' % control_cond_name)
    C = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % test_cond_name)
    D = clustio.ParseNormal('auc_results/%s_results_reweight_RAW.txt' % control_cond_name)
    assert len(C) == len(D)
    for i, j in comb(xrange(len(C)), 2):
        vC = N.sqrt(C.M[i][j] * C.M[j][i])
        vD = N.sqrt(D.M[i][j] * D.M[j][i])
        C.M[i][j] = C.M[j][i] = vC
        D.M[i][j] = D.M[j][i] = vD
    A.add_method('M8002', 'Enhanced in %s vs Normal' % test_cond_name, 'C', color='purple')
    A.add_method('M8003', 'Enhanced in Normal vs %s' % test_cond_name, 'C', color='green')
    A.add_method('M7000', 'Abnormal Link', 'C', color='red')
    A.add_method('M7001', 'Deleted Link', 'C', color='lightblue')
    Anodes = A.nodes
    And = dict([ (x.name, x) for x in Anodes ])
    comp = clustio.ParseNormal('auc_results/%s_vs_%s_thresh_95.txt' % (test_cond_name, control_cond_name))
    nodenames = And.keys()
    for i, j in comb(xrange(len(nodenames)), 2):
        n1 = nodenames[i]
        n2 = nodenames[j]
        ilen = len(pathway_dict[n1.lower()])
        jlen = len(pathway_dict[n2.lower()])
        thresh = comp.get_samples([str(ilen)]).get_features([str(jlen)]).M[0][0]
        wB = B.get_samples([n1.lower()]).get_features([n2.lower()]).M[0][0]
        wA = 0.0
        if A.isconnected(n1, n2):
            wA = float([ edge.weight for edge in And[n1].links if edge.target == n2 ][0])
        wC = C.get_samples([n1.lower()]).get_features([n2.lower()]).M[0][0]
        wD = D.get_samples([n1.lower()]).get_features([n2.lower()]).M[0][0]
        meth = None
        if wA and wB:
            if thresh and (wC - wD) >= thresh:
                meth = 'M8002'
            elif thresh and (wD - wC) >= thresh:
                meth = 'M8003'
            else:
                meth = 'M0099'
        elif wA and not wB:
            if thresh and (wC - wD) >= thresh:
                meth = 'M7000'
            else:
                meth = 'M0099'
        elif wB and not wA:
            meth = 'M7001'
        if meth is not None:
            if meth == 'M7001':
                A.add_edge(And[n1], And[n2], method=meth)
            else:
                link1 = [ link for link in And[n1].links if link.target == n2 ][0]
                link2 = [ link for link in And[n2].links if link.target == n1 ][0]
                link1.method = meth
                link2.method = meth
    return A

def prettify_tree(tree):

    for node in tree.nodes:
        node.ncc = 'white'
        node.labelStyle = '1'

        # KEGG pathway style
        node.data.type = '6'
        node.childVisible = 'false'
        # We might consider labeling them for the user as well

        for link in node.links:
            link.weight = '0.2'

    return tree

def create_tree(test_cond_name, control_cond_name, pathway_dict, infile_suffix=NO_COMPARISON_SUFFIX, outfile_suffix=COMPARISON_SUFFIX):

    print
    print('Creating VisML tree %s_%s.xml...' % (test_cond_name, infile_suffix))
    tree = scripts.create_visml_from_sdata(clustio.ParseNormal('sig_connections/%s_sig_connections_95.txt' % test_cond_name), dmatrix=True)
    tree = prettify_tree(tree)
    tree.write('%s_%s.xml' % (test_cond_name, infile_suffix))

    print('Creating VisML tree %s_%s.xml...' % (test_cond_name, outfile_suffix))
    tree = difftree(pathway_dict, test_cond_name, control_cond_name, infile_suffix=infile_suffix)
    tree.write('%s_%s.xml' % (test_cond_name, outfile_suffix))
