#!/usr/bin/python
"""Dispatch pairwise dependency computations between two matrices.

USE EXAMPLE:
python dispatch_pairwise.py outdir=/fs/lustre/osu6683/gse15745 tabfile_1=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab tabfile_2=$HOME/gse15745/GSE15745.GPL8178.miRNA.normed.tab dry=True
"""
from util import *
import sys
import shutil



BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s start=%(start)d end=%(end)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

def dispatch_pairwise(tabfile_1=None, tabfile_2=None, outdir=WORK_DIR, function=None, k=200000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=2, n_ppn=12, walltime='6:00:00'):
  assert tabfile_1, tabfile_2
  n_nodes, n_ppn, start_offset, k = map(int, (n_nodes, n_ppn, start_offset, k))
  assert k > 1 and start_offset >= 0 and n_nodes > 0 and n_ppn > 0
  if jobname is None: 
    jobname = os.path.basename(tabfile)
  
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

  npy_fname_1, n1 = npy_varlist_from_tabfile(tabfile_1, outdir)
  npy_fname_2, n2 = npy_varlist_from_tabfile(tabfile_2, outdir)

  # Move npy_fname to work_dir
  work_npy_fname_1 = move_numpy_to_workdir(work_dir, npy_fname_1)
  work_npy_fname_2 = move_numpy_to_workdir(work_dir, npy_fname_2)

  # dispatch jobs in a loop
  num_pairs = int(n * (n-1) / 2) # no diagonal: n choose 2

  # Write jobs to dispatch script in a list.
  t = tstamp()
  for function in all_functions:
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
        'stdout_fname': os.path.join(work_dir, "log_%s_%s_%s.out" % (jobname, t, function)),
        'stderr_fname': os.path.join(work_dir, "log_%s_%s_%s.err" % (jobname, t, function)),
        'function': function,
      }
      fp.write(cmd); fp.write('\n')
      offset += k
    fp.close()
    # Submit job script.
    qsub_script = QSUB_TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, 'dispatch_script': dispatch_script_fname}
    print qsub_script
    if not dry:
      fork_qsub_submit(qsub_script)
      print "Batch job submitted."
    else:
      print "Dry run."


if __name__ == "__main__":
  dispatch_pairwise(**dict([s.split('=') for s in sys.argv[1:]]))
