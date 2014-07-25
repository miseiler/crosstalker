# File output utilities for ConsensusCluster

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

def write_normal(sdata, filename):
    """

    Takes an sdata obj and writes out a tab-delimited datafile, suitable for ParseNormal
    useful to convert between data formats, writing PCA-selected data, etc

    """
    
    if not len(sdata.gene_names) > 0: 
        raise ValueError, "No gene names found! Unsuitable for this data format."

    sids = [x.sample_id for x in sdata.samples]

    f = open(filename, 'w')

    #Sample line, first row
    f.write("\t".join(['SAMPLE ID'] + sids))
    f.write("\n")

    #Data
    for i in xrange(len(sdata.gene_names)):
        f.write("\t".join([sdata.gene_names[i]] + [ str(sdata.M[j][i]) for j in xrange(len(sdata.samples)) ]))
        f.write("\n")

    f.close()

def write_table(ndict, filename):
    """Write a tab delimited flat file, one key per line"""

    ls = ndict.keys()
    ls.sort()

    f = open(filename, 'w')
    
    for key in ls:
        f.write("\t".join([str(key), str(ndict[key])]))
        f.write("\n")

    f.close()

                        
def write_list(ls, filename):
    """

    Write a simple list to a file, each element on one line.

    Suitable for cluster definitions, etc

    """

    f = open(filename, 'w')

    for x in ls:
        f.write(str(x))
        f.write('\n')

    f.close()

