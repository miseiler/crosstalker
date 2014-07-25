"""

Text parsers for creating SampleData objects


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

import sys
import numpy as N
import ioutils

def read_cluster_log(log):
    """Get clusters in dict format from log"""

    clust_data = {}
    current_cluster = None
    
    f = open(log, 'r')

    for line in f:
        spl = line.strip().split()

        if spl and spl[0] == 'Cluster':
            current_cluster = spl[1] + " " + spl[2]
            clust_data[current_cluster] = []

        elif line and line[0] == '\t' and current_cluster is not None:
            tspl = line.strip().split('\t')
            clust_data[current_cluster].append(tspl[0])

    return clust_data

def read_table(filename):
    """

    Read a simple conversion table between one name and another, useful for converting between probe names and gene symbols.

    This really should parse the entire CSV in the future, but memory concerns have held me back for now.
    Maybe an SQLite database?

    Returns a dict of first-column: second-column associations

    """

    conv = dict()

    handle = open(filename, 'r')

    for line in handle:
        a = line.split("\t")

        if a:

            conv[a[0]] = " - ".join(a[1:]).strip()

    handle.close()

    return conv

def read_list(filename):
    """

    read_list(filename)

        Simply returns a concatenated list of every line in filename, with the newlines stripped

    Returns: A list of strings

    """

    handle = open(filename, 'r')
    entries = [ x.strip() for x in handle ]
    handle.close()

    return entries


class SampleData(object):
    """

    BaseCluster

        Usage:
            
            BaseParserObj.samples.append(SampleData(sample_id, sample_num, sample_class, data, index))
        
            sample_id       - Label for the sample.  Need not be unique. (Optional, but highly recommended)
            sample_num      - Another label (Optional)
            sample_class    - Yet another label, usually for known subclasses (Optional)
            index           - Yet another label (Optional)

        Properties

            cluster_id      - Id of cluster to which this sample belongs.  Generally assigned by clustering algorithms.


    """

    def __init__(self, sample_id=None, sample_num=None, sample_class=None, cluster_id=None, index=None):

        self.cluster_id = None
        self.sample_class = sample_class
        self.sample_id = sample_id
        self.sample_num = sample_num
        self.index = index

    def __getitem__(self, x):

        return self.sample_id[x]


class BaseParser(object):
    """

    BaseParser

        Common methods and actions for all parsers
        Designed to be subclassed.  All parsers must do the following:
        
        1: Have an __init__ method which takes data_file as an argument

        2: Have a _parse_data_file method which adds SampleData objects to self.samples

        3: If the data vector entries have labels (e.g., gene probe names), they should be assigned
           to self.gene_names as a list, in the same length and order as the data vectors in each sample
    
        4: Be called "ParseSOMETHING" where SOMETHING is any name
           The gtk front end looks for ParseX names and puts them in a drop-down list for convenience

        Properties

            samples     - A list of SampleData objects
            gene_names  - A list of data vector entry labels

    """

    def __init__(self, data_file):

        self.samples = []  #List of SampleData instances
        self.gene_names = []

        self._parse_data_file(data_file)

    def _parse_data_file(self, data_file):
        """Parse datafile into sample name<->number pairs and load data"""
        pass

    def __delitem__(self, x):

        inds = list(xrange(len(self.samples)))
        inds.pop(x)

        self.M = self.M.take(tuple(inds), 0)
        
        del self.samples[x]

    def __getitem__(self, x):

        return self.samples[x]

    def __len__(self):

        return len(self.samples)

    def __iter__(self):

        return iter(self.samples)

    def __add__(self, x):
        """
        Concatenate two BaseParser classes. This will concatenate samples if the gene sets are identical, and vice versa.
        Currently requires that the sample/gene names be in the same order
        This will try to self.merge if these conditions are not satisfied.
        
        """

        c = NullParser()

        if len(self) == len(x) and self.sample_ids == x.sample_ids:

            c.samples = self.samples
            c.M = N.hstack([self.M, x.M])
            c.gene_names = N.hstack([self.gene_names, x.gene_names])
            return c

        elif len(self.gene_names) == len(x.gene_names) and (self.gene_names == x.gene_names).all():

            x_sids = set(x.sample_ids)

            for sid in self.sample_ids:
                if sid in x_sids:
                    raise ValueError, 'One or more sample ids are not unique, unable to concatenate'

            c.samples = self.samples + x.samples
            c.M = N.vstack([self.M, x.M])
            c.gene_names = self.gene_names
            return c

        return self.merge(x)

    def get_samples(self, *filenames):
        """
        
        Extract a new BaseParser object from samples in *filenames

        Filenames can be a list, a dict, or filenames.

        If a dict or filenames, get_samples will intelligently set the sample_class attribute
        on each sample to the dict key or file name pertaining to it. See ioutils.list_or_files
        for detail.

        The new BaseParser object will be reduced to the sample set in *filenames.

        """

        sample_defs    = ioutils.list_or_files(*filenames)

        # The list of samples is a special case, and we shouldn't mangle sample_classes in that case
        if sample_defs.keys() == [None]:
            d = self.class_dict
        else:
            d = sample_defs
        
        class_dict = dict([ (d[k][i], k) for k in d
                            for i in xrange(len(d[k])) ])

        sample_indices = tuple(ioutils.argintersect(self.sample_ids, sum(sample_defs.values(), [])))

        c = NullParser()

        c.samples = [ SampleData(cluster_id=self.samples[i].cluster_id, sample_id=self.samples[i].sample_id,
                                 sample_num=self.samples[i].sample_num, index=self.samples[i].index,
                                 sample_class=class_dict[self.samples[i].sample_id]) for i in sample_indices ]

        c.M = self.M.take(sample_indices, 0)
        c.gene_names = self.gene_names.copy()

        return c

    def get_features(self, *filenames):
        """

        Extract a new BaseParser object from features in *filenames

        Filenames can be a list, a dict, or filenames. The complete feature set will be concatenated
        from all available keys or file contents. See ioutils.list_or_files for detail.

        The new BaseParser object will be reduced to the feature set in *filenames.

        """

        plist = sum(ioutils.list_or_files(*filenames).values(), [])
        pinds = tuple(ioutils.argintersect(self.gene_names, plist))

        c = NullParser()

        c.M = self.M.take(pinds, 1)
        c.gene_names = self.gene_names.take(pinds)
        c.samples = [ SampleData(cluster_id=x.cluster_id, sample_id=x.sample_id,
                                 sample_num=x.sample_num, index=x.index,
                                 sample_class=x.sample_class) for x in self.samples ]

        return c

    def get_class(self, cls):
        """Extract a new BaseParser object composed of all samples with the sample_class cls"""

        # We can't use argintersect because of the uniqueness constraint
        if cls in self.sample_classes:

            sample_indices = tuple([ i for i in xrange(len(self)) if self.samples[i].sample_class == cls ])

            c = NullParser()

            c.samples = [ SampleData(cluster_id=self.samples[i].cluster_id, sample_id=self.samples[i].sample_id,
                                     sample_num=self.samples[i].sample_num, index=self.samples[i].index,
                                     sample_class=cls) for i in sample_indices ]

            c.M = self.M.take(sample_indices, 0)
            c.gene_names = self.gene_names.copy()

            return c

        raise ValueError, '%s not found in dataset' % (cls)

    def merge(self, x):
        """Merge two BaseParser objects. The feature set is reduced to the most common set, and samples are appended."""

        x_sids = dict.fromkeys(x.sample_ids)

        for sid in self.sample_ids:
            if sid in x_sids:
                raise ValueError, 'One or more sample ids are not unique, unable to concatenate'

        if len(self.gene_names) and len(x.gene_names):
            fst = ioutils.argintersect(self.gene_names, x.gene_names)
            snd = ioutils.argintersect(x.gene_names, self.gene_names)

            if not len(fst):
                raise ValueError, "No matching gene names in either set! Cannot concatenate."

            gene_names = self.gene_names.take(tuple(fst))
            
            M = self.M.take(fst, 1) #Memory-intensive but we certainly don't want to
            Q = x.M.take(snd, 1)    #change the objects themselves

        else:
            print('WARNING: Concatenating one or more sets without a gene list! Do so at your own risk!')

        c = NullParser()

        c.gene_names = gene_names
        c.samples = self.samples + x.samples
        c.M = N.vstack((M, Q))

        return c

    @property
    def sample_ids(self):
        """Returns a list of sample_ids from self.samples"""

        return [ x.sample_id for x in self.samples ]

    @property
    def sample_classes(self):
        """Returns a sorted list of sample_classes from self.samples"""

        return list(set([ x.sample_class for x in self.samples ]))
    
    @property
    def class_dict(self):
        """Returns a dict with sample_classes as keys and the representative sample_ids as values"""

        cd = {}
        for sam in self.samples:
            cd.setdefault(sam.sample_class, []).append(sam.sample_id)
        return cd

class NullParser(BaseParser):
    """

    Null class.  Used to create new BaseParser objects.

    """

    def __init__(self):
        pass


class ParseNormal(BaseParser):
    """

    ParseNormal

        This is the one you use when the data is just a table with no special characteristics
        Probes in rows, samples in columns.  Sample ids are in row 0, and gene ids
        are in column 0.
        
    """

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_file):
        """Parse datafile into sample name<->number pairs and load probe data"""

        sample_list = open(data_file, 'r').readline().strip('\n').split('\t')

        for sam_id in sample_list[1:]:
            self.samples.append(SampleData(sample_id=sam_id))

        self.gene_names = N.loadtxt(data_file, delimiter='\t', skiprows=1, usecols=[0], dtype='S')

        self.M = N.loadtxt(data_file, delimiter='\t', skiprows=1, usecols=xrange(1, len(sample_list)), dtype=N.float32).T

        # Fixes 1-sample bug
        if len(self.M.shape) == 1:
            self.M = self.M.reshape(1, self.M.shape[0])

        try:
            assert len(set(self.sample_ids)) == len(self)
        except:
            raise ValueError, 'One or more sample ids in this file are not unique!'

        try:
            assert len(set(self.gene_names)) == len(self.gene_names)
        except:
            raise ValueError, 'One or more features in this file are not unique!'
