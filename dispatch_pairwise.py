#!/usr/bin/python
"""Dispatch pairwise dependency computations.

Note: .npy now should be .pkl

USE EXAMPLE:
python dispatch_pairwise.py outdir=/fs/lustre/osu6683/gse15745 tabfile=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab dry=True

time python $HOME/dependency_matrix_dispatcher/dispatch_pairwise.py outdir=/fs/lustre/osu6683/gse15745 npy_file=/nfs/01/osu6683/gse15745/old/GSE15745.GPL6104.mRNA.normed.tab.npy_aligned_with_GSE15745.GPL8490.methyl.normed.tab.npy.npy.TCTX.npy dry=True
"""
from util import *
import sys
import shutil
import math
import cPickle as pickle


BATCH_CMD = "time python %(script_path)s/%(script_name)s npyfile=%(npyfile)s start=%(start)d end=%(end)d work_dir=%(work_dir)s function=%(function)s n=%(n)d > %(stdout_fname)s 2> %(stderr_fname)s"

def dispatch_pairwise(tabfile=None, pkl_file=None, outdir=WORK_DIR, function=None, k=200000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=6, n_ppn=12, walltime='24:00:00'):
  assert bool(tabfile is None) != bool(pkl_file is None)
  if pkl_file:
    assert pkl_file.rpartition('.')[2].lower() == "pkl"
  n_nodes, n_ppn, start_offset, k = map(int, (n_nodes, n_ppn, start_offset, k))
  assert k > 1 and start_offset >= 0 and n_nodes > 0 and n_ppn > 0
  if jobname is None:
    if tabfile is not None:
      jobname_base = os.path.basename(tabfile)
    else:
      jobname_base = os.path.basename(pkl_file)
  else:
    jobname_base = jobname
  
  if function is None:
    all_functions = FUNCTIONS
  else:
    assert function in FUNCTIONS
    all_functions = {function: FUNCTIONS[function]}

  # create outdir if it does not exist
  if not os.path.exists(outdir):
    print "Make outdir %s..." % (outdir)
    make_dir(outdir)
  else:
    print "outdir %s exists." % (outdir)
  # chdir to outdir
  os.chdir(outdir)
  print "Changed working directory to outdir %s." % (outdir)

  # for each function, create subdirs
  outdirs = {}
  for function in all_functions:
    path = os.path.join(outdir, function)
    make_dir(path)
    outdirs[function] = os.path.abspath(path)

  # If provided a pickled numpy matrix already, then we don't need to covert it.
  if pkl_file:
    print "Given npy matrix %s. Do not convert from .tab text file." % (pkl_file)
    M = pickle.load(open(pkl_file))
    n = np.size(M,0)
    print "M is (%d by %d)" % (np.size(M,0), np.size(M,1))
    del M
  # Convert .tab into npy matrix
  else:
    npy_fname, n, M = npy_varlist_from_tabfile(tabfile, outdir)
    del M

  # Move npy_fname to workdir
  bool_copy = (tabfile is None) # if not provided tabfile, then copy npy, don't move it
  if bool_copy:
    print "Tab file not converted; copy npy matrix rather than move it."
  work_npy_fname = move_numpy_to_workdir(work_dir, pkl_file, bool_copy)
  
  # dispatch jobs in a loop
  num_pairs = int(n * (n-1) / 2) # no diagonal: n choose 2

  # Write jobs to dispatch script in a list.
  t = tstamp()
  for function in all_functions:
    # Handle custom batch scripts
    if type(all_functions[function]) == str:
      print "%s type is string '%s', not function. Assume that it is the name of custom batch file." \
          % (function, all_functions[function])
      script_fname = all_functions[function]
    else:
      script_fname = "batch_pairwise.py"
      
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
        'script_fname': script_fname,
        'npyfile': work_npy_fname,
        'start': start,
        'end': end,
        'n': n,
        'work_dir': outdirs[function],
        'stdout_fname': os.path.join(outdirs[function], "log_%s_%s_%s_%d_%d.out" % (jobname, t, function, start, end)),
        'stderr_fname': os.path.join(outdirs[function], "log_%s_%s_%s_%d_%d.err" % (jobname, t, function, start, end)),
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
