#!/usr/bin/python
import subprocess
import datetime
from dcor import *
import os
import errno
import numpy.ma as ma
from scipy.stats import mstats
import shutil

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
mpiexec parallel-command-processor %(dispatch_script)s
"""

def euclidean(x,y):
  q=x-y
  return ma.sqrt((q*q.T).sum())

FUNCTIONS = {
  'pearson': lambda x, y: mstats.pearsonr(x,y)[0],
  'spearman': lambda x, y: mstats.spearmanr(x,y)[0],
  'euclidean': euclidean,
  'kendalltau': lambda x,y: mstats.kendalltau(x,y)[0],
  'dcor': dcor,
  }


def move_numpy_to_workdir(work_dir, npy_fname):
  work_npy_fname = os.path.join(work_dir, os.path.basename(npy_fname))
  if not os.path.exists(work_npy_fname):
    shutil.move(npy_fname, work_npy_fname)
    print "Moved %s to %s." % (npy_fname, work_npy_fname)
  return work_npy_fname

def read_samples(tabfile_1_coltitles):
  sample_titles_1 = [s for s in open(tabfile_1_coltitles).next().strip('\n').split('\t')]
  sample_titles_set_1 = set(sample_titles_1)
  idx_1 = dict([(i, s) for i, s in enumerate(sample_titles_1)])
  return sample_titles_1, sample_titles_set_1, idx_1


def npy_varlist_from_tabfile(tabfile, outdir):
  """Import tabfile into numpy matrix and variable list; save to disk.

  Args:
    tabfile: str of path to .tab input data file
  Returns:
    (str, int) of path to saved numpy object and number of rows in data matrix
  """
  varlist_fname = os.path.join(outdir, os.path.basename(tabfile) + ".varlist.txt")
  npy_fname = os.path.join(outdir, "%s.npy" % tabfile)
  if os.path.exists(varlist_fname) and os.path.exists(npy_fname):
    print "Both %s and %s exist, do not recreate varlist and numpy masked matrix files." % \
        (varlist_fname, npy_fname)
    # load numpy matrix to get its size
    M = ma.load(npy_fname)
    n = np.size(M, 0)
  else:
    # import tab
    varlist = []
    # M is masked matrix
    print "Loading %s into masked numpy matrix and varlist..." % tabfile
    M = np.genfromtxt(name_iter(open(tabfile), varlist), usemask=True, delimiter='\t')
    n = np.size(M, 0) # number of rows (variables)
    # save to file
    print "Saving matrix and varlist..."
    fp = open(varlist_fname, "w")
    fp.write('\n'.join(varlist))
    fp.close()
    ma.dump(M, npy_fname)
  return npy_fname, n, M
  


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
    return True

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

    
