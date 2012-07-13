#!/usr/bin/python
from util import *
import numpy.ma as ma
from py_symmetric_matrix import *

#BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 1000

def main(npyfile=None, offset=None, work_dir=None, function=None, n=None, start=None, end=None, batchname=None):
  assert npyfile, work_dir
  assert function in FUNCTIONS
  n = int(n)
  assert n > 0 and k > 0 and offset >= 0
  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npyfile), function, start, end)

  if start is None:
    start = 0
  else:
    start = int(start)
  if end is None: 
    end = n*(n-1) / 2
  else:
    end = int(end)

  M = ma.load(npyfile)
  f = FUNCTIONS[function]
  R = np.zeros(end-start)
  print "Starting to write %d pairs for %s" % (end-start, batchname)
  for i, j in enumerate(xrange(start, end)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(i, m)
    R[i] = f(M[x], M[y])
  print "Computed %d pairs for %s" % (end-start, batchname)

  output_fname = os.path.join(work_dir, batchname+".txt")
  print "Saving results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, R)
  print "Saved %s." % output_fname
  
  
