#!/usr/bin/python
import subprocess
import datetime
from dcor import *
import os
import errno
import numpy.ma as ma
from scipy.stats import mstats


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

    
