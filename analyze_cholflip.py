#!/usr/bin/python
# Code by Peter Kasson, 2015

"""Analyze sterol flips across a bilayer."""

import commands
import gflags
import numpy
import os
import sys

import clusteratoms
from HierTools import xvg_parse


# Overall idea is take a phosphate index spec and a sterol oxygen index spec
# Bicluster the phosphates -> Z positions
# Use clusteratoms.py for this
# then assign sterol oxygen to nearest
# monitor changes

def assign_chol(pdbfile, trajfile, chol_ndx, leaflet_ndx, gmx_path='/usr/local/gromacs/bin'):
  """Assigns sterol atoms to leaflets in each frame of trajectory.
  Args:
    pdbfile:  pdb, gro, or tpr file
    trajfile:  trajectory file
    chol_ndx:  index file with one group:  sterol atoms to assign
    leaflet_ndx:  index file with one group per leafleta
    gmx_path:  path to Gromacs binaries
  Rets:
    chol_assigns:  array with rows = timepoints, cols = sterols, val = leaflet.
  """
  # make joint ndx file
  tmp_ndx = '%s.ndx' % str(uuid.uuid4())
  os.system('cat %s %s > %s' % (chol_ndx, leaflet_ndx, tmp_ndx))
  # count number of groups in leaflet index
  nleaflets = int(commands.getoutput('grep -c "\[" %s' % leaflet_ndx))
  tmp_xvg = '%s.xvg' % str(uuid.uuid4())
  os.system('%s/g_mindist -s %s -f %s -n %s -or %s -ng %d -respertime'
            % (gmx_path, pdbfile, trajfile, tmp_ndx, tmp_xvg, nleaflets))
  os.unlink('mindist.xvg')
  # read xvg back in
  # first N cols should be first leaflet distance, then second leaflet
  leaflet_dist = xvg_parse(tmp_xvg, None, keycol=0)
  # assign
  chol_assigns = numpy.reshape(leaflet_dist, (nrows, nleaflets, -1)).argmin(1)
  os.unlink(tmp_ndx)
  os.unlink(tmp_xvg)
  return chol_assigns


def analyze_flips(cholmat):
  """Given a matrix of leaflet assignments, find flips.
  Args:
    cholmat:  array with rows = timepoints, cols = sterols, val = which leaflet
  Rets:
    flipcount:  array with number of flips per sterol.
  """
  return [sum(numpy.equal(x, numpy.roll(x, 1))) for x in cholmat.T]


def chol_flips(pdbfile, trajfile, lipid_atomname, chol_selstring,
                   gmx_path = '/usr/local/gromacs/bin'):
  """Make index files, run analysis on cholesterol flips.
  Args:
    pdbfile:  pdb file
    trajfile:  trajectory file
    lipid_atomname:  regexp for lipid selection (leaflet clustering)
    chol_selstring:  selection string for sterol atoms (one per residue)
    gmx_path:  path to Gromacs binaries
  Rets:
  """
