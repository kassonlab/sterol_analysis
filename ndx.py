#
#  ndx.py
#
#
#  Written by Peter Kasson, 2011-2012
#  Copyright (c) 2012

"""Tools for manipulating index files."""

def write_index(index_dict, fname):
  """Writes a Gromacs index file.
  Args:
    index_dict:  dict describing index
    fname:  filename
  """
  #based on Gromacs code.  converts 0-indexed to 1-indexed
  outfile = open(fname, 'w')
  ngroups = len(index_dict['IndexNames'])
  for groupidx in range(ngroups):
    outfile.write('[ %s ]\n' % index_dict['IndexNames'][groupidx])
    linectr = 0
    for val in index_dict['IndexGroups'][groupidx]:
      outfile.write('%4d ' % (val+1))
      if linectr % 15 == 14:
        outfile.write('\n')
      linectr += 1
    outfile.write('\n')
  outfile.close()
