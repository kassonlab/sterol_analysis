#!/usr/bin/python
# Code by Peter Kasson, 2014

"""Agglomerative clustering of atoms in a structure."""

import re
import sys

import gflags
import numpy
from scipy.cluster import hierarchy

import ndx
import PDB


def clusteratoms(PDBfilename, atomnamesel, nclusters, writendx=None, resname='',
                 radial_dist=False):
  """Cluster atoms into k groups by nearest neighbor.
    Args:
      PDBfilename:  input file
      atomnamesel: regular expression for atom names to select
      nclusters: number of clusters to make
      writendx: optional filename of Gromacs index to write
    Rets:
      clusterdata: array where each row is atomid, resid, clusternum.
  """
  infile = open(PDBfilename, 'r')
  pdbrecord = PDB.readPDB(infile)
  atomlines = []
  idxctr = 0
  for line in pdbrecord[0]:
    if line.__class__ in [PDB.ATOM, PDB.HETATM]:
      if re.match(atomnamesel, line.name) and re.match(resname, line.resName):
        atomlines.append([line.serial, line.resSeq, idxctr,
                          line.x, line.y, line.z])
    idxctr += 1
  atomarr = numpy.array(atomlines)
  print atomarr.shape
  # alternate approach:
  # specify maxdist, then:
  # clusteridx = hierarchy.fclusterdata(atomarr[:, 3:6], maxdist,
  #                                     criterion='distance')
  if radial_dist:
    ctr_coord = numpy.mean(atomarr[:, 3:6], 0)
    radial_arr = numpy.array([[numpy.linalg.norm(ctr_coord - line[3:6])]
                              for line in atomarr])
    clusteridx = hierarchy.fclusterdata(radial_arr, nclusters,
                                        criterion='maxclust')
  else:
    clusteridx = hierarchy.fclusterdata(atomarr[:, 3:6], nclusters,
                                        criterion='maxclust')
  # slightly kludgy construction of return data structure
  clusterdata = numpy.zeros([len(clusteridx), 3])
  clusterdata[:, 0] = atomarr[:, 0]
  clusterdata[:, 1] = atomarr[:, 1]
  clusterdata[:, 2] = clusteridx
  # option to write a Gromacs index
  if writendx:
    indexdict = {}
    indexdict['IndexNames'] = ['Cluster%d' % i for i in range(nclusters)]
    indexdict['IndexGroups'] = []
    for i in range(nclusters):
      indexdict['IndexGroups'].append(atomarr[numpy.nonzero(clusteridx == i+1)[0], 2])
    ndx.write_index(indexdict, writendx)
  return clusterdata


if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('infile', '',
                       'Input PDB')
  gflags.DEFINE_string('atomsel', '(.+)',
                       'Regexp to select atom names')
  gflags.DEFINE_string('ressel', '',
                       'Optional regex to select atom names')
  gflags.DEFINE_string('outfile', 'out.dat',
                       'Output file name')
  gflags.DEFINE_string('ndxfile', None,
                       'Optional output index')
  gflags.DEFINE_integer('numclusters', 2,
                       'Number of clusters to make')
  argv = FLAGS(sys.argv)
  clusters = clusteratoms(FLAGS.infile, FLAGS.atomsel, FLAGS.numclusters,
                          FLAGS.ndxfile)
  numpy.savetxt(FLAGS.outfile, clusters, fmt='%d\t%d\t%d')
