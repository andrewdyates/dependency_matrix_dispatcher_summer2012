#!/usr/bin/python
"""Compile a directory of blocks pairwise numpy matrices into a single squareform matrix.

GSE15745.GPL6104.mRNA.normed.tab.npy_dcor_0_200000.npy

SAMPLE USE:

python $HOME/dependency_matrix_dispatcher/compile_pairwise.py path=/fs/lustre/osu6683/gse15745/dcor outpath_prefix=$HOME/gse15745/gse15745_gpl6104_dcor n=246053836

TODO: automatically compute n
"""
from __future__ import division
from py_symmetric_matrix import *
import re
import numpy as np
import os, sys

RX = re.compile('.*?_(\d+)_(\d+)\.npy')


def main(path, outpath_prefix, n=None, npy_fname=None):
  assert path, outpath_prefix

  if n is not None:
    n = int(n)
  else:
    assert npy_fname is not None
    M = np.load(npy_fname)
    n = np.size(M,0)
  assert n > 0

  M = np.zeros(n, dtype=np.float32)
  B = np.zeros(n, dtype=np.bool)
  n_set_total, n_dupe_total, n_nan_total = 0, 0, 0
  for fname in os.listdir(path):
    m = RX.match(fname)
    if m:
      start, end = int(m.group(1)), int(m.group(2))
      Q = np.load(os.path.join(path,fname))
      n_set, n_dupe, n_nan = 0, 0, 0
      for i, x in enumerate(range(start, end)):
        M[x] = Q[i]
        if np.isnan(Q[i]):
          n_nan += 1
        elif not B[x]:
          B[x] = 1
          n_set += 1
        else:
          n_dupe += 1
      print "%.2f%% Complete: Set %d (%d dupes, %d nan) from %s. Expected %d." % \
          ((end-start)/n_set*100, n_set, n_dupe, n_nan, fname, end-start)
      n_set_total += n_set
      n_dupe_total += n_dupe
      n_nan_total += n_nan
  print "%.2f%% Complete. Set %d (%d dupes, %d nan) from %s. Expected %d." % (n_set_total/n*100, n_set_total, n_dupe_total, n_nan_total, path, n)
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
