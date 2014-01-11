#!/usr/bin/python
# Code by Peter Kasson, 2014

"""Agglomerative clustering of atoms in a structure."""

import re
import sys

import gflags
import numpy
from scipy.cluster import hierarchy

import PDB


def clusteratoms(PDBfilename, atomnamesel, nclusters):
  """Cluster atoms into k groups by nearest neighbor.
    Args:
      PDBfilename:  input file
      atomnamesel: regular expression for atom names to select
      nclusters: number of clusters to make
    Rets:
      clusterdata: array where each row is atomid, resid, clusternum.
  """
  infile = open(PDBfilename, 'r')
  pdbrecord = PDB.readPDB(infile)
  atomlines = []
  for line in pdbrecord[0]:
    if line.__class__ in [PDB.ATOM, PDB.HETATM]:
      if re.search(atomnamesel, line.name):
        atomlines.append([line.serial, line.resSeq, line.x, line.y, line.z])
  atomarr = numpy.array(atomlines)
  clusteridx = hierarchy.fclusterdata(atomlines[:, 2:4], nclusters,
                                      criterion='maxclust')
  clusterdata = atomarr[:, 0:2]
  clusterdata[:, 2] = clusteridx
  return clusterdata


if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('infile', '',
                       'Input PDB')
  gflags.DEFINE_string('atomsel', '(.+)',
                       'Regexp to select atom names')
  gflags.DEFINE_string('outfile', 'out.dat',
                       'Output file name')
  gflags.DEFINE_integer('numclusters', 2,
                       'Number of clusters to make')
  argv = FLAGS(sys.argv)
  clusters = clusteratoms(FLAGS.infile, FLAGS.atomsel, FLAGS.numclusters)
  numpy.savetxt(FLAGS.outfile, clusters)
