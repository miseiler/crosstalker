# Common IO Utilities for ConsensusCluster

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
import parsers

def new_filename(filename, default, extension):
    """Try to avoid squashing old files by incrementing -num extensions to the filename"""

    if filename is None:
        filename = default

    newname = filename + extension
    dirlist = os.listdir(os.curdir)
    counter = 0

    while newname in dirlist:
        counter += 1
        newname = "".join([filename, ' - %s' % counter, extension])

    return newname

def list_or_files(*args):
    """
    
    Return the contents of *args as dict[name] = list pairs regardless of whether the user provided a filename or dict
    If a list is provided, it will be appended to the dict with key ''

    Returns a dict of name->list pairs

    """

    ndict = {}

    for lst in args:
        if isinstance(lst, list):
            ndict.setdefault(None, []).extend(lst)
            #print("WARNING: List received! Multiple lists are concatenated!")

        elif isinstance(lst, dict):
            ndict.update(lst)
        
        else:
            name = os.path.basename(lst)
            ndict.setdefault(name, []).extend(parsers.read_list(lst))

    return ndict

def union_old(list1, list2):
    """
    
    DEPRECATED
    
    Return the indices which make up the union between two lists in O(n) time.

    Returns a tuple, where tuple[0] is the list of indices in list1 which is in common with list2, and tuple[1]
    is the same list for list2

    Assumes that both lists are composed of unique items.

    """

    try:
        assert len(set(list1)) == len(list1) and len(set(list2)) == len(list2)
    except:
        raise ValueError, 'One or more sets contain non-unique elements'

    swapped = False
    if len(list1) > len(list2):         #Make list2 the longer one
        list1, list2 = list2, list1
        swapped = True

    indices_list1 = N.argsort(list1)
    indices_list2 = N.argsort(list2)

    union_indices_list1 = []
    union_indices_list2 = []
    
    breakpoint = 0

    for i in indices_list1:    
        for j in range(len(indices_list2))[breakpoint:]:    #Ugly, but reduces complexity
            idx = indices_list2[j]

            if list1[i] == list2[idx]:
                union_indices_list1.append(i)
                union_indices_list2.append(idx)
                breakpoint = j
                break

    if not swapped:
        return union_indices_list1, union_indices_list2

    return union_indices_list2, union_indices_list1

def argintersect(list1, list2):
    """

    Returns the indices of list1 members which exist in list2

    Throws an error if the lists are not unique

    """

    try:
        assert len(set(list1)) == len(list1) and len(set(list2)) == len(list2)
    except:
        raise ValueError, 'One or more sets contain non-unique elements'

    # Oct 23 2012
    # Bug in numpy 1.7 (perhaps others?) which in very specific cases prevents
    # in1d from returning anything if searched data structure is not pre-cast to array
    list1 = N.array(list1)

    return [ x[0] for x in N.argwhere(N.in1d(list1, list2, assume_unique=True)) ]


def get_indices(s, filename):
    """
    Return the indices of the samples in filename in the sdata object
    
    """

    sams = list_or_files(filename)
    name = sams.keys()[0]
    
    return name, argintersect([ x.sample_id for x in s.samples ], sams[name])

