#!/usr/bin/python
"""Compute dependencies between two matrices.

1) load both matrices
2) get offset row of first matrix
3) compute all dependencies of offset row with other matrix
4) save results to disk

Only compute dCOR for pairs with both values present.
"""
from util import *
import numpy.ma as ma
import sys

#BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 1000

def main(npyfile_1=None, npyfile_2=None, offset=None, work_dir=None, function=None, batchname=None):
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
  f = FUNCTIONS[function]
  if function == "dcor":
    print "dCOR implementation:", DCOR_LIB
  n = np.size(M2, 0)
  R = np.zeros(n)
  n_nan = 0
  print "Starting to write %d pairs for %s" % (n, batchname)
  for i in xrange(n):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, n, batchname)
    # Mask pairs with at least one missing value
    shared_mask = ~(M1[offset].mask | M2[i].mask)
    R[i] = f(M1[offset][shared_mask], M2[i][shared_mask])

    if np.isnan(R[i]):
      n_nan += 1
      m1_masked = np.count_nonzero(M1.mask)
      m2_masked = np.count_nonzero(M2.mask)
      print "NAN!", offset, i, m1_masked, m2_masked
    
  print "Computed %d pairs for %s using %s." % (n, batchname, function)
  print "%d nans" % (n_nan)

  output_fname = os.path.join(work_dir, batchname+".npy")
  print "Saving %d results at %s." % (n, output_fname)
  np.save(output_fname, R)
  print "Saved %s." % output_fname
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
