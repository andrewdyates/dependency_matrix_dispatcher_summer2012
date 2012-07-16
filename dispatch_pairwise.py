#!/usr/bin/python
"""Dispatch pairwise dependency computations.

USE EXAMPLE:
python dispatch_pairwise.py outdir=/fs/lustre/osu6683/gse15745 tabfile=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab dry=True
"""
from util import *
import sys
import shutil
import math


BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s start=%(start)d end=%(end)d work_dir=%(work_dir)s function=%(function)s n=%(n)d > %(stdout_fname)s 2> %(stderr_fname)s"

def dispatch_pairwise(tabfile=None, outdir=WORK_DIR, function=None, k=200000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=6, n_ppn=12, walltime='8:00:00'):
  assert tabfile
  n_nodes, n_ppn, start_offset, k = map(int, (n_nodes, n_ppn, start_offset, k))
  assert k > 1 and start_offset >= 0 and n_nodes > 0 and n_ppn > 0
  if jobname is None: 
    jobname_base = os.path.basename(tabfile)
  else:
    jobname_base = jobname
  
  if function is None:
    all_functions = FUNCTIONS
  else:
    assert function in FUNCTIONS
    all_functions = {function: FUNCTIONS[function]}

  # for each function, create subdirs
  outdirs = {}
  for function in all_functions:
    path = os.path.join(outdir, function)
    make_dir(path)
    outdirs[function] = os.path.abspath(path)

  # Convert .tab into npy matrix
  npy_fname, n, M = npy_varlist_from_tabfile(tabfile, outdir)
  del M

  # Move npy_fname to workdir
  work_npy_fname = move_numpy_to_workdir(work_dir, npy_fname)
  
  # dispatch jobs in a loop
  num_pairs = int(n * (n-1) / 2) # no diagonal: n choose 2

  # Write jobs to dispatch script in a list.
  t = tstamp()
  for function in all_functions:
    jobname = "%s_%s" % (jobname_base, function)
    # Create new dispatch script per function
    dispatch_script_fname = \
      make_script_name(work_dir, os.path.basename(work_npy_fname), "dispatch_%s" % function)
    print "Creating batch script '%s'..." % dispatch_script_fname
    fp = open(dispatch_script_fname, 'w')
    offset = start_offset
    while offset < num_pairs:
      start = offset
      end = offset + k
      if end > num_pairs:
        end = num_pairs
      cmd = BATCH_CMD % {
        'script_path': os.path.dirname(os.path.realpath(__file__)),
        'npyfile': work_npy_fname,
        'start': start,
        'end': end,
        'n': n,
        'work_dir': outdirs[function],
        'stdout_fname': os.path.join(work_dir, "log_%s_%s_%s_%d_%d.out" % (jobname, t, function, start, end)),
        'stderr_fname': os.path.join(work_dir, "log_%s_%s_%s_%d_%d.err" % (jobname, t, function, start, end)),
        'function': function,
      }
      fp.write(cmd); fp.write('\n')
      offset += k
    fp.close()
    num_jobs = math.ceil(num_pairs / k)
    # Submit job script.
    qsub_script = QSUB_TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, \
      'dispatch_script': dispatch_script_fname, 'num_pairs': num_pairs, 'k': k, 'num_jobs': num_jobs, 'function': function}
    print qsub_script
    if not dry:
      fork_qsub_submit(qsub_script)
      print "Batch job submitted."
    else:
      print "Dry run."


if __name__ == "__main__":
  dispatch_pairwise(**dict([s.split('=') for s in sys.argv[1:]]))
