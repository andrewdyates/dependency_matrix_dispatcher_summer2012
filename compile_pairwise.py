#!/usr/bin/python
"""Compile a directory of blocks pairwise numpy matrices into a single squareform matrix.

Files in same path are expected to be in format:

[id]_[start]_[end].[matrix].npy

where:
  id: identifying name of source data matrix
  start: int of start index
  end: int of end index+1
  matrix: name of statistic saved in matrix
  
examples:
  GSE7307_GPL570.normed.masked.pkl.gt0.25.pina.tab.pkl_dcor_5220000_5250000.DCOR.npy
  GSE7307_GPL570.normed.masked.pkl.gt0.25.pina.tab.pkl_dcor_5220000_5250000.DCOV.npy

SAMPLE USE:
  python $HOME/dependency_matrix_dispatcher/compile_pairwise.py path=/fs/lustre/osu6683/gse15745/dcor outpath_prefix=$HOME/gse15745/gse15745_gpl6104_dcor n=246053836
"""
from __future__ import division
from py_symmetric_matrix import *
import re
import numpy as np
import os, sys

RX = re.compile('(?P<id>.*?)_(?P<start>\d+)_(?P<end>\d+)\.(?P<matrix>\w+).npy')

class Result(object):
  def __init__(self, n, dtype=np.float):
    assert n > 0
    self.M = np.zeros(n, dtype=dtype)
    self.B = np.zeros(n, dtype=np.bool)
    self.n_set_total, self.n_dupe_total, self.n_nan_total = 0, 0, 0
    self.Q_last = None

def main(path, outpath_prefix, n=None, npy_fname=None, precision=32):
  assert path, outpath_prefix
  precision = int(precision)
  assert precision in (32, 64)

  if n is not None:
    n = int(n)
  else:
    assert npy_fname is not None
    M = np.load(npy_fname)
    n = np.size(M,0)
  assert n > 0

  if precision == 32:
    print "Compiling to precision np.float32."
    dtype=np.float32
  elif precision == 64:
    print "Compiling to precision np.float64 (np.float)."
    dtype=np.float64
  else:
    print "Unknown precision:", precision
    sys.exit(1)

  # We don't know in advance how many matrices we will need to fill.
  Results = {}
  for fname in os.listdir(path):
    m = RX.match(fname)
    if m:
      # Parse file name into parts.
      start = int(m.group('start'))
      end = int(m.group('end'))
      matrix_name = m.group('matrix')
      Q = np.load(os.path.join(path, fname))
      
      # Get Result for this matrix_name
      if matrix_name not in Results:
        print "Creating new Result matrix '%s' of size %d, type %s." % (matrix_name, n, str(dtype))
        Results[matrix_name] = Result(n=n, dtype=dtype)
      R = Results[matrix_name]
        
      # Check to make sure that Q is not the same as last Q for this Result matrix.
      if R.Q_last is not None:
        try:
          assert np.size(Q) != np.size(R.Q_last) or np.abs(np.sum(Q - R.Q_last)) > 0
        except AssertionError:
          print "Matrix segment seems repeated..."
          print fname
          print Q[:5]
          print R.Q_last[:5]
          print (Q-R.Q_last)[:5]
          print np.sum(Q-R.Q_last)
          raise
      n_set, n_dupe, n_nan = 0, 0, 0
      for i, x in enumerate(range(start, end)):
        R.M[x] = Q[i]
        if np.isnan(Q[i]):
          # Implicit: B[x] = 0
          n_nan += 1
        elif not R.B[x]:
          R.B[x] = 1
          n_set += 1
        else:
          n_dupe += 1
      R.Q_last = Q
      print "%.2f%% Complete: Set %d (%d dupes, %d nan) from %s. Expected %d." % \
          ((end-start)/n_set*100, n_set, n_dupe, n_nan, fname, end-start)
      R.n_set_total += n_set
      R.n_dupe_total += n_dupe
      R.n_nan_total += n_nan

  # Save results for each result.
  for matrix_name, R in Results.items():
    print "Saving results for matrix %s..." % (matrix_name)
    print "%.2f%% Complete. Set %d (%d dupes, %d nan, %d unmasked) from %s. Expected %d." % \
        (R.n_set_total/n*100, R.n_set_total, R.n_dupe_total, R.n_nan_total, np.sum(R.B), path, n)
    M_fname = "%s.%s.values.npy" % (outpath_prefix, matrix_name)
    B_fname = "%s.%s.isset.npy" % (outpath_prefix, matrix_name)
    
    np.save(M_fname, R.M)
    print "Saved %s." % (M_fname)
    if R.n_set_total != n:
      np.save(B_fname, R.B)
      print "!!! Because values are missing, saved %s." % (B_fname)
    else:
      print "No values missing; did not save boolean 'isset' matrix."

  print "Compilation of %s complete. Saved %d result matrices." % (path, len(Results))

  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
