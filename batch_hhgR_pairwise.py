#!/usr/bin/python
"""Special batch script for HHG in R.
Adapted from `batch_pairwise.py`

Normalizes scores from 0 to 1 by maximum score; see HHG_R package readme.

Computes (and saves matrices for):
  * sum_chisquared
  * sum_logratio
  * max_chisquared
  * max_logratio

Unlike other batch_pairwise.py script, this computes several matrices at once; one for
  each of the HHG statistics.

EXAMPLE MANUAL USE:
  time python $HOME/dependency_matrix_dispatcher/batch_hhgR_pairwise.py npyfile=/fs/lustre/osu6683/GSE7307.normed.tab.pkl work_dir=/fs/lustre/osu6683/gse7307/hhgR_test n=54675 start=0 end=10 batchname=hhgR_test verbose=True
"""
from util import *
import numpy as np
import cPickle as pickle
import numpy.ma as ma
from py_symmetric_matrix import *
import sys
# R and Rpy2 must be installed.
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
REPORT_N = 50000

def main(npyfile=None, work_dir=None, n=None, start=None, end=None, batchname=None, verbose=False, *args, **kwds):
  function = "hhgR"
  assert npyfile, work_dir
  n, start, end = map(int, (n, start, end))
  assert n > 0 and start >= 0 and end > 0
  if type(verbose) == str and verbose.lower() in ("none", "f", "false"):
    verbose = False
  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npyfile), function, start, end)
  # Load HHG library from R installation.
  r('library("HHG2x2")')

  if verbose:
    live_test()

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
  # load level 2 pickle rather than numpy's inefficient MaskedArray ASCII binary format if .pkl
  if npyfile.rpartition('.')[2].lower() == 'pkl':
    print "Loading as level 2 pickle"
    M = pickle.load(open(npyfile))
  else:
    print "Loading as level 0 numpy.MaskedArray pickle"
    M = ma.load(npyfile)
  SUM_CHI, SUM_LR, MAX_CHI, MAX_LR = np.zeros(end-start), np.zeros(end-start), np.zeros(end-start), np.zeros(end-start)
  n_nan_sum_chi, n_nan_sum_lr, n_nan_max_chi, n_nan_max_lr = 0, 0, 0, 0
  print "Starting to write %d pairs for %s" % (end-start, batchname)
  for i, j in enumerate(xrange(start, end)):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(i, n)
    assert x >= 0 and y >= 0
    # Mask values with at least one missing value in pair.
    try:
      shared_mask = ~(M[x].mask | M[y].mask)
    except IndexError:
      print "WARNING! INDEX ERROR! i %d, x %d, y %d, n %d" %(i,x,y,n)
      raise
    n_values = np.sum(shared_mask)

    # Put vector pair in R namespace
    robjects.globalenv["x"] = M[x][shared_mask].data
    robjects.globalenv["y"] = M[y][shared_mask].data
    # Execute HHG algorithm in R.
    r('Dx = as.matrix(dist((x),diag=TRUE,upper=TRUE))')
    r('Dy = as.matrix(dist((y),diag=TRUE,upper=TRUE))')
    HHG = r('myHHG(Dx,Dy)')
    # Extract results from R namespace; normalize
    SUM_CHI[i], SUM_LR[i], MAX_CHI[i], MAX_LR[i] = \
        float(HHG.rx('sum_chisquared')[0][0]) / ((n_values)*(n_values-2)*(n_values-3)), \
        float(HHG.rx('sum_lr')[0][0]) / ((n_values)*(n_values-2)*(n_values-3)), \
        float(HHG.rx('max_chisquared')[0][0]) / (n_values-2), \
        float(HHG.rx('max_lr')[0][0]) / (n_values-2)/np.log(2)
    if np.isnan(SUM_CHI[i]):
      n_nan_sum_chi += 1
    if np.isnan(SUM_LR[i]):
      n_nan_sum_lr += 1
    if np.isnan(MAX_CHI[i]):
      n_nan_max_chi += 1
    if np.isnan(MAX_LR[i]):
      n_nan_max_lr += 1
    if verbose:
      print "%d: %d %.4f %.4f %.4f %.4f" % (i, n_values, SUM_CHI[i], SUM_LR[i], MAX_CHI[i], MAX_LR[i])
      
  print "Computed %d pairs for %s" % (end-start, batchname)
  print "%d sum_chi nans, %d sum_lr nans, %d max_chi nans, %d max_lr nans" % \
      (n_nan_sum_chi, n_nan_sum_lr, n_nan_max_chi, n_nan_max_lr)
  n_bad = sum((n_nan_sum_chi, n_nan_sum_lr, n_nan_max_chi, n_nan_max_lr))
  if n_bad > 0:
    print "!!!WARNING: There exists at least one (%d) not-a-numbers (nans) in this batch." % n_bad

  # Save each of 4 matrices
  output_fname = os.path.join(work_dir, batchname+".sum_chi.npy")
  print "Saving SUM CHI results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, SUM_CHI)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".sum_lr.npy")
  print "Saving SUM LOG RATIO results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, SUM_LR)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".max_chi.npy")
  print "Saving MAX CHI results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MAX_CHI)
  print "Saved %s." % output_fname

  output_fname = os.path.join(work_dir, batchname+".max_lr.npy")
  print "Saving MAX LOG RATIO results %d through %d as %s. (zero-indexed)" % (start, end-1, output_fname)
  np.save(output_fname, MAX_LR)
  print "Saved %s." % output_fname


def live_test(n_values=80):
  print "Test run for identity function n=%d..." % (n_values)
  robjects.globalenv["x"] = np.arange(n_values)
  robjects.globalenv["y"] = np.arange(n_values)
  r('Dx = as.matrix(dist((x),diag=TRUE,upper=TRUE))')
  r('Dy = as.matrix(dist((y),diag=TRUE,upper=TRUE))')
  HHG = r('myHHG(Dx,Dy)')
  print HHG
  SUM_CHI_TEST, SUM_LR_TEST, MAX_CHI_TEST, MAX_LR_TEST = \
      float(HHG.rx('sum_chisquared')[0][0]) / ((n_values)*(n_values-2)*(n_values-3)), \
      float(HHG.rx('sum_lr')[0][0]) / ((n_values)*(n_values-2)*(n_values-3)), \
      float(HHG.rx('max_chisquared')[0][0]) / (n_values-2), \
      float(HHG.rx('max_lr')[0][0]) / (n_values-2)/np.log(2)
  print "Max at %d: %.4f %.4f %.4f %.4f" % (n_values, SUM_CHI_TEST, SUM_LR_TEST, MAX_CHI_TEST, MAX_LR_TEST)
  
  
if __name__ == "__main__":
  print sys.argv
  main(**dict([s.split('=') for s in sys.argv[1:]]))
