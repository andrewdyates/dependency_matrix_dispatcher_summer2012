#!/usr/bin/python
from util import *
import numpy.ma as ma
from py_symmetric_matrix import *
import sys
import cPickle as pickle

#BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 50000

def main(npyfile=None, work_dir=None, function=None, n=None, start=None, end=None, batchname=None):
  
  assert npyfile, work_dir
  assert function in FUNCTIONS
  n, start, end = map(int, (n, start, end))
  assert n > 0 and start >= 0 and end > 0
  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npyfile), function, start, end)

  # Do not recreate existing files.
  output_fname = os.path.join(work_dir, batchname+".npy")
  if os.path.exists(output_fname):
    print "%s already exists. Exiting..." % output_fname
    return 1

  if start is None:
    start = 0
  else:
    start = int(start)
  if end is None: 
    end = n*(n-1) / 2
  else:
    end = int(end)

  print "Loading %s..." % npyfile
  if npyfile.rpartition('.')[2].lower() == 'pkl':
    print "Loading as level 2 pickle"
    M = pickle.load(open(npyfile))
  else:
    print "Loading as level 0 numpy.MaskedArray pickle"
    M = ma.load(npyfile)
  M = ma.load(npyfile)
  f = FUNCTIONS[function]
  if function == "dcor":
    print "dCOR implementation:", DCOR_LIB
  R = np.zeros(end-start)
  n_nan = 0
  print "Starting to write %d pairs for %s" % (end-start, batchname)
  for i, j in enumerate(xrange(start, end)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(i, n)
    assert x >= 0 and y >= 0
    # TODO: mask missing values like shared_mask = ~(M1[offset].mask | M2[i].mask)
    try:
      shared_mask = ~(M[x].mask | M[y].mask)
      R[i] = f(M[x][shared_mask], M[y][shared_mask])
    except IndexError:
      print "WARNING! INDEX ERROR! i %d, x %d, y %d, n %d" %(i,x,y,n)
      raise
    if np.isnan(R[i]):
      n_nan += 1

  print "Computed %d pairs for %s" % (end-start, batchname)
  print "%d nans" % (n_nan)

  output_fname = os.path.join(work_dir, batchname+".npy")
  print "Saving results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, R)
  print "Saved %s." % output_fname
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
