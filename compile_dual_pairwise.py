#!/usr/bin/python
"""Compile a directory of dual pairwise numpy matrix rows into a single matrix.

SAMPLE USE:

python $HOME/dependency_matrix_dispatcher/compile_dual_pairwise.py path=/fs/lustre/osu6683/gse15745_gpl8178_gpl6104/dcor outpath_prefix=$HOME/gse15745_gpl8178_gpl6104_dcor n_rows=735 n_cols=22184
"""
from __future__ import division
from py_symmetric_matrix import *
import re
import numpy as np
import os, sys

RX = re.compile('.*?_(\d+)\.npy')


def main(path, outpath_prefix, n_rows, n_cols):
  assert path, outpath_prefix
  n_rows, n_cols = int(n_rows), int(n_cols)
  assert n_rows > 0 and n_cols > 0
  n = n_cols*n_rows

  M = np.zeros((n_rows, n_cols), dtype=np.float32)
  B = np.zeros((n_rows, n_cols), dtype=np.bool)
  
  n_set_total, n_dupe_total, n_nan_total = 0, 0, 0
  for fname in os.listdir(path):
    m = RX.match(fname)
    if m:
      row_num = int(m.group(1))
      Q = np.load(os.path.join(path,fname))
      q_nans = np.count_nonzero(np.isnan(Q))
      assert q_nans == 0, q_nans
      n_set, n_dupe, n_nan = 0, 0, 0
      assert np.size(Q,0) == n_cols
      for i in xrange(n_cols):
        M[row_num,i] = Q[i]
        if np.isnan(Q[i]):
          n_nan += 1
        elif not B[row_num,i]:
          B[row_num,i] = 1
          n_set += 1
        else:
          n_dupe += 1
        # assert not np.isnan(M[row_num,i])
      print "%.2f%% Complete: Set %d (%d dupes, %d n_nan) from %s. Expected %d." % \
          ((n_set)/n_cols*100, n_set, n_dupe, n_nan, fname, n_cols)
      n_set_total += n_set
      n_dupe_total += n_dupe
      n_nan_total += n_dupe
      
  print "%.2f%% Complete. Set %d (%d dupes, %d nan) from %s. Expected %d." % (n_set_total/n*100, n_set_total, n_dupe_total, n_nan_total, path, n)
  assert np.count_nonzero(np.isnan(M)) == 0
  
  M_fname, B_fname = outpath_prefix+".values.npy", outpath_prefix+".isset.npy"
  np.save(M_fname, M)
  print "Saved %s." % (M_fname)
  if n_set_total != n:
    np.save(B_fname, B)
    print "Because values are missing, saved %s." % (B_fname)
  else:
    print "No values missing; did not save boolean 'isset' matrix."

  print "Compilation of %s complete." % path

  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
