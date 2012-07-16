#!/usr/bin/python
"""Compute dependencies between two matrices.

1) load both matrices
2) get offset row of first matrix
3) compute all dependencies of offset row with other matrix
4) save results to disk
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

  M1 = ma.load(npyfile_1)
  M2 = ma.load(npyfile_2)
  assert offset < np.size(M1, 0)
  f = FUNCTIONS[function]
  n = np.size(M2, 0)
  R = np.zeros(n)
  print "Starting to write %d pairs for %s" % (n, batchname)
  for i in xrange(n):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, n, batchname)
    try:
      R[i] = f(M[x], M[y])
    except IndexError:
      print "WARNING! INDEX ERROR! i %d, x %d, y %d, n %d" %(i,x,y,n)
      raise
  print "Computed %d pairs for %s using %s." % (n, batchname, function)

  output_fname = os.path.join(work_dir, batchname+".npy")
  print "Saving %d results at %s." % (n, output_fname)
  np.save(output_fname, R)
  print "Saved %s." % output_fname
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
