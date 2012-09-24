#!/usr/bin/python
"""Compute dependencies between two matrices.
Compare row "offset" of matrix M1 with all other rows of matrix M2.

TODO: this should be integrated with batch_pairwise to reduce redundant code.

1) load both matrices
2) get offset row of first matrix
3) compute all dependencies of offset row with other matrix
4) save results to disk
"""
from util import *
import numpy.ma as ma
import sys

#BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 50000

def main(npyfile_1=None, npyfile_2=None, offset=None, work_dir=None, function=None, batchname=None, verbose=False):
  assert npyfile_1 and npyfile_2 and work_dir
  assert function in FUNCTIONS
  offset = int(offset)
  assert offset >= 0
  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_vs_%s_%s_%d" % \
      (os.path.basename(npyfile_1), os.path.basename(npyfile_2), function, offset)

  M1 = ma.load(npyfile_1) # these should be pickles
  M2 = ma.load(npyfile_2)
  assert offset < np.size(M1, 0)
  assert np.count_nonzero(np.isnan(M1.compressed())) == 0
  assert np.count_nonzero(np.isnan(M2.compressed())) == 0

  # Get batch function handler for this function.
  size = np.size(M2,0)
  F = FUNCTIONS[function](size)
  if verbose:
    print "Vebose: Try F on n=700 identity."
    print F.compute_one(np.arange(700), np.arange(700))

  # Compute pairs using batch handler `F`
  print "Starting to write %d pairs for %s" % (size, batchname)
  for i in xrange(np.size(M2,0)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, size, batchname)
    # Mask pairs with at least one missing value
    shared_mask = ~(M1[offset].mask | M2[i].mask)
    X, Y = M1[offset][shared_mask].data, M2[i][shared_mask].data
    assert np.size(X) == np.size(Y) <= np.size(M1,1)

    F.compute(X,Y,i)
    if verbose:
      d = F.get(i)
      s = " ".join(["%s=%f"%(k,v) for k,v in d.items()])
      print "%d: " % (i), s
    
  print "Computed %d pairs for %s using %s." % (n, batchname, function)
  n_nans = F.nans()
  print "%d nans" % (n_nans)
  if n_nans > 0:
    print "!!!WARNING: There exists at least one (%d) not-a-numbers (nans) in this batch." % (n_nans)

  out_names = F.save(work_dir, batchname)
  print "Saving %d results as:" % (n)
  for name, out_name in out_names.items():
    print "%s: %s" % (name, out_name)
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
