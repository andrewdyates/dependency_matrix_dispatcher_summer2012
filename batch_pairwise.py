#!/usr/bin/python
"""Compute batch of pairwise dependencies.

EXAMPLE MANUAL USE:
  time python $HOME/dependency_matrix_dispatcher/batch_pairwise.py function=dcor npyfile=/fs/lustre/osu6683/GSE7307.normed.tab.pkl work_dir=/fs/lustre/osu6683/gse7307/manual_test n=54675 start=5 end=15 verbose=True
"""
from util import *
import numpy.ma as ma
from py_symmetric_matrix import *
import sys
import cPickle as pickle

#BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 50000

def main(npyfile=None, work_dir=None, function=None, n=None, start=None, end=None, batchname=None, verbose=False, *args, **kwds):
  
  assert npyfile and work_dir
  assert function in FUNCTIONS
  if start is None:
    start = 0
  else:
    start = int(start)
  if end is None: 
    end = n*(n-1) / 2
  else:
    end = int(end)
  n = int(n)
  assert n > 0 and start >= 0 and end > 0

  # If batchname not provided, compute default value. Use batchname in output file name.
  if batchname is None or batchname.lower() in ("false", "f", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npyfile), function, start, end)

  # Do not recreate existing batch output files.
  output_fname = os.path.join(work_dir, batchname+".npy")
  if os.path.exists(output_fname):
    print "%s already exists. Exiting..." % output_fname
    return 1
  
  # Load data file
  print "Loading %s..." % npyfile
  if npyfile.rpartition('.')[2].lower() == 'pkl':
    print "Loading as level 2 pickle"
    M = pickle.load(open(npyfile))
  else:
    print "Loading as level 0 numpy.MaskedArray pickle"
    M = ma.load(npyfile)

  # TODO: Revise this as some kind of factory selector.
  f = FUNCTIONS[function]
  if function == "dcor":
    print "dCOR implementation:", DCOR_LIB

  # Create vector containers
  R = np.zeros(end-start)
  print "Starting to write %d pairs for %s" % (end-start, batchname)
  for i, j in enumerate(xrange(start, end)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(j, n)
    assert x >= 0 and y >= 0 and x < np.size(M,0) and y < np.size(M,0)
    idx_check = sym_idx(x,y,n)
    assert idx_check == j and idx_check >= start and idx_check < end

    shared_mask = ~(M[x].mask | M[y].mask)
    X, Y = M[x][shared_mask].data, M[y][shared_mask].data
    assert np.size(X) == np.size(Y) <= np.size(M,1)

    
    R[i] = f(X,Y)
    if verbose:
      print i, j, x, y, R[i]

  n_nan = np.sum(np.isnan(R))
  print "Computed %d pairs for %s" % (end-start, batchname)
  print "%d nans" % (n_nan)
  if n_nan > 0:
    print "!!!WARNING: There exists at least one (%d) not-a-numbers (nans) in this batch." % (n_nan)

  output_fname = os.path.join(work_dir, batchname+".npy")
  print "Saving results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, R)
  print "Saved %s." % output_fname
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
