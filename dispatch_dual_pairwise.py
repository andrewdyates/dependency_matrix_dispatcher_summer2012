#!/usr/bin/python
"""Dispatch pairwise dependency computations between two matrices.

USE EXAMPLE:
python dispatch_pairwise.py outdir=/fs/lustre/osu6683/gse15745 tabfile_1=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab tabfile_2=$HOME/gse15745/GSE15745.GPL8178.miRNA.normed.tab dry=True
"""
from util import *
import sys
import shutil


BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"

def dispatch_pairwise(tabfile_1=None, tabfile_2=None, tabfile_1_coltitles=None, tabfile_2_coltitles=None, outdir=WORK_DIR, function=None, k=200000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=2, n_ppn=12, walltime='6:00:00'):
  assert tabfile_1 and tabfile_2 and tabfile_1_coltitles and tabfile_2_coltitles
  n_nodes, n_ppn, start_offset, k = map(int, (n_nodes, n_ppn, start_offset, k))
  assert k > 1 and start_offset >= 0 and n_nodes > 0 and n_ppn > 0
  if jobname is None: 
    jobname = os.path.basename("%s_vs_%s" % (tabfile_1, tabfile_2))
  
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

  npy_fname_1, n1, M1 = npy_varlist_from_tabfile(tabfile_1, outdir)
  npy_fname_2, n2, M2 = npy_varlist_from_tabfile(tabfile_2, outdir)
  # let the smaller matrix be matrix 1
  if n1 > n2:
    npy_fname_1, npy_fname_2 = npy_fname_2, npy_fname_1
    tabfile_1_coltitles, tabfile_2_coltitles = tabfile_2_coltitles, tabfile_1_coltitles
    n1, n2 = n2, n1
    M1, M2 = M2, M1
    print "Swapped matrices 1 and 2 so that matrix 1 is smaller."

  # Column-Align two matrices.
  titles1, set1, idx1 = read_samples(tabfile_1_coltitles)
  titles2, set2, idx2 = read_samples(tabfile_2_coltitles)
  cols = []
  for s in titles1:
    if s in set2:
      cols.append((idx1[s], idx2[s]))
  print "Aligned %d cols from %s and %d cols from %s to %d shared columns." % \
    (len(titles1), tabfile_1_coltitles, len(titles2), tabfile_2_coltitles, len(cols))

  # TODO:
  # reorder matrices and write them to file
    
  # Move npy_fname to work_dir
  work_npy_fname_1 = move_numpy_to_workdir(work_dir, npy_fname_1)
  work_npy_fname_2 = move_numpy_to_workdir(work_dir, npy_fname_2)

  # Write jobs to dispatch script in a list.
  t = tstamp()
  for function in all_functions:
    # Create new dispatch script per function
    dispatch_script_fname = \
      make_script_name(work_dir, os.path.basename(work_npy_fname_1), "dispatch_%s" % function)
    print "Creating batch script '%s'..." % dispatch_script_fname
    fp = open(dispatch_script_fname, 'w')
    offset = start_offset
    for offset in xrange(n1):
# BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"
      cmd = BATCH_CMD % {
        'script_path': os.path.dirname(os.path.realpath(__file__)),
        'npyfile_1': work_npy_fname_1,
        'npyfile_2': work_npy_fname_2,
        'offset': offset,
        'work_dir': outdirs[function],
        'stdout_fname': os.path.join(work_dir, "log_%s_%s_%s.out" % (jobname, t, function)),
        'stderr_fname': os.path.join(work_dir, "log_%s_%s_%s.err" % (jobname, t, function)),
        'function': function,
      }
      fp.write(cmd); fp.write('\n')
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