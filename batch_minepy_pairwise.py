#!/usr/bin/python
"""Special batch script for cmine at http://mpba.fbk.eu/cmine.
Adapted from `batch_pairwise.py`

Computes (and saves matrices for):
  mic, mas, mev, mcn

Unlike other batch_pairwise.py script, this computes several matrices at once; one for
  each of the MINE statistics.
"""
from util import *
import numpy as np
import numpy.ma as ma
from py_symmetric_matrix import *
import sys
# minepy must be installed. See http://mpba.fbk.eu/cmine
import minepy

#BATCH_CMD = "time python %(script_path)s/batch_minepy_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

REPORT_N = 50000

def main(npyfile=None, work_dir=None, n=None, start=None, end=None, batchname=None, alpha=0.6, c=15):
  function = "minepy"
  assert npyfile, work_dir
  n, start, end = map(int, (n, start, end))
  assert n > 0 and start >= 0 and end > 0
  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npyfile), function, start, end)
  alpha = float(alpha)
  c = int(c)
  assert alpha > 0 and alpha < 1
  assert c >= 0

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

  print "Loading %s..." % (npyfile)
  # TODO: load level 2 pickle rather than numpy's inefficient MaskedArray ASCII binary format.
  M = ma.load(npyfile)
  MIC, MAS, MEV, MCN = np.zeros(end-start), np.zeros(end-start), np.zeros(end-start), np.zeros(end-start)
  n_nan_mic, n_nan_mas, n_nan_mev, n_nan_mcv = 0, 0, 0, 0
  print "Starting to write %d pairs for %s" % (end-start, batchname)
  for i, j in enumerate(xrange(start, end)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(i, n)
    assert x >= 0 and y >= 0
    # Create minepy computation object
    mine = minepy.MINE(alpha=alpha, c=c)
    try:
      mine.score(M[x], M[y])
    except IndexError:
      print "WARNING! INDEX ERROR! i %d, x %d, y %d, n %d" %(i,x,y,n)
      raise
    MIC[i], MAS[i], MEV[i], MCV[i] = mine.mic(), mine.mas(), mine.mev(), mine.mcv()
    if np.isnan(MIC[i]):
      n_nan_mic += 1
    if np.isnan(MAS[i]):
      n_nan_mas += 1
    if np.isnan(MEV[i]):
      n_nan_mev += 1
    if np.isnan(MCV[i]):
      n_nan_mcv += 1

  print "Computed %d pairs for %s" % (end-start, batchname)
  print "%d mic nans, %d mas nans, %d mev nans, %d mcv nans" % (n_nan_mic, n_nan_mas, n_nan_mev, n_nan_mcv)
  n_bad = sum((n_nan_mic, n_nan_mas, n_nan_mev, n_nan_mcv))
  if n_bad > 0
    print "!!!WARNING: There exists at least one (%d) not-a-numbers (nans) in this batch." % n_bad

  # Save each of 4 matrices
  output_fname = os.path.join(work_dir, batchname+".mic.npy")
  print "Saving MIC results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MIC)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".mas.npy")
  print "Saving MAS results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MAS)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".mev.npy")
  print "Saving MEV results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MEV)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".mcv.npy")
  print "Saving MCV results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MCV)
  print "Saved %s." % output_fname

  
if __name__ == "__main__":
  print sys.argv
  main(**dict([s.split('=') for s in sys.argv[1:]]))
