#!/usr/bin/python
"""Note: this file contains redundant functionality implemented in other modules."""
import subprocess
import datetime
import cPickle as pickle
import os
import errno
from numpy import ma
from scipy.stats import mstats
import shutil

MISSING_VALUES = ",".join(["nan", "None", "N/A", "none", "NaN"])

WORK_DIR = "/fs/lustre/osu6683"
QSUB_TEMPLATE = \
"""#PBS -N %(jobname)s
#PBS -l nodes=%(n_nodes)d:ppn=%(n_ppn)d
#PBS -j oe
#PBS -S /bin/bash
#PBS -l walltime=%(walltime)s
#tdate=$(date +%%T)

set -x
cd /nfs/01/osu6683/
source .bash_profile
echo "Dispatching %(num_jobs)d jobs each with %(k)d pairs from %(num_pairs)s total pairs."
echo "Function: %(function)s"
mpiexec parallel-command-processor %(dispatch_script)s
"""

from batch import *
from dcor_batch import *
from mine_batch import *
from hhg_batch import *


FUNCTIONS = {
  'pearson': PCCBatch,
  'covariance': CovBatch,
  'spearman': SpearmanBatch,
  'dcor': DcorBatch,
  'mine': MINEBatch,
  'hhg': HHGBatch,
  'euclidean': EuclideanBatch,
  'kendalltau': KendallBatch,
  }
IGNORE = ['euclidean', 'kendalltau']

def move_numpy_to_workdir(work_dir, npy_fname, do_copy=False):
  work_npy_fname = os.path.join(work_dir, os.path.basename(npy_fname))
  if not os.path.exists(work_npy_fname):
    if do_copy:
      shutil.copyfile(npy_fname, work_npy_fname)
      print "Copied %s to %s." % (npy_fname, work_npy_fname)
    else:
      shutil.move(npy_fname, work_npy_fname)
      print "Moved %s to %s." % (npy_fname, work_npy_fname)
  return work_npy_fname

def read_samples(tabfile_1_coltitles):
  sample_titles_1 = [s for s in open(tabfile_1_coltitles).next().strip('\n').split(',')]
  assert len(sample_titles_1) > 1, "delimit sample titles with commas"
  sample_titles_set_1 = set(sample_titles_1)
  idx_1 = dict([(s, i) for i, s in enumerate(sample_titles_1)])
  return sample_titles_1, sample_titles_set_1, idx_1


def npy_varlist_from_tabfile(tabfile, outdir, overwrite=False):
  """Convert .tab into pickled npy np.ma.MaskedArray matrix and variable list; save to disk.

  Args:
    tabfile: str of path to .tab input data file
  Returns:
    (str, int, np.ma.MaskedArray) of path to pickled numpy object, rows in data matrix, M
  """
  varlist_fname = os.path.join(outdir, os.path.basename(tabfile) + ".varlist.txt")
  npy_fname = os.path.join(outdir, "%s.pkl" % tabfile)
  if not overwrite and os.path.exists(varlist_fname) and os.path.exists(npy_fname):
    print "Both %s and %s exist, do not recreate varlist and pickled numpy masked matrix files." % \
        (varlist_fname, npy_fname)
    # load numpy matrix to get its size
    try:
      M = pickle.load(open(npy_fname))
    except EOFError:
      # Some error may have occurred in a previous save. Reload .tab and overwrite pickle.
      print "Error reading %s. Attempt to recreate and overwrite." % (npy_fname)
      return npy_varlist_from_tabfile(tabfile, outdir, overwrite=True)
    n = np.size(M, 0)
  else:
    # import tab
    varlist = []
    # M is masked matrix
    print "Loading %s into masked numpy matrix and varlist..." % tabfile
    M = np.genfromtxt(name_iter(open(tabfile), varlist), usemask=True, delimiter='\t', missing_values=MISSING_VALUES, dtype=np.float)
    n = np.size(M, 0) # number of rows (variables)
    # save to file
    print "Saving %d rows in varlist to %s..." % (len(varlist), varlist_fname)
    fp = open(varlist_fname, "w")
    fp.write('\n'.join(varlist))
    fp.close()
    print "Saving %d rows, %d columns of matrix as protocol 2 pickle to %s..." % \
        (np.size(M,0), np.size(M,1), npy_fname)
    pickle.dump(M, open(npy_fname, 'w'), protocol=2)
    
  n_nans = np.count_nonzero(np.isnan(M.compressed()))
  assert n_nans == 0, "%d 'nan's exists in matrix %s!" % (n_nans, npy_fname)
  return npy_fname, n, M
  

def print_matrix_stats(M1):
  n_missing1 = np.count_nonzero(M1.mask)
  t1 = np.size(M1.ravel())
  print "Missing values: %d of %d (%d by %d) (%.2f%%)" % \
      (n_missing1, t1, np.size(M1,0), np.size(M1,1), n_missing1/t1)
  print "nans: Total: %d, With Mask: %d" % (np.count_nonzero(np.isnan(M1)), np.count_nonzero(np.isnan(M1.compressed())))

def make_dir(outdir):
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

def tstamp():
  return datetime.datetime.isoformat(datetime.datetime.now())
    
def make_script_name(work_dir, tab_basename, process_name):
  tmp_script_name = "tmp_script_%s_%s_%s.sh" % (tab_basename, process_name, tstamp())
  dispatch_script_fname = os.path.join(work_dir, tmp_script_name)
  return dispatch_script_fname

def is_dry(dry):
  if type(dry) == str and dry.lower() in ('false', 'f', 'none', ''):
    return False
  else:
    return bool(dry)

def count_tab_rows(tabfile):
  x = len([None for f in open(tabfile) if f[0] not in ('\n', '#')])
  return x

def fork_qsub_submit(qsub_script):
  p = subprocess.Popen("qsub", stdin=subprocess.PIPE)
  p.communicate(input=qsub_script)
  p.stdin.close()

def name_iter(fp, varlist):
  """Load a labeled row matrix as a line iterator."""
  for line in fp:
    if line[0] in ('#', '\n'): continue
    name,c,row = line.partition('\t')
    varlist.append(name)
    yield row

    
