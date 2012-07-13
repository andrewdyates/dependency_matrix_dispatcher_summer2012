#!/usr/bin/python
from util import *
import sys
import shutil


BATCH_CMD = "time python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s n=%(n)d >> %(stdout_fname)s 2>> %(stderr_fname)s"

def dispatch_pairwise(tabfile=None, outdir=None, function=None, k=500000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=2, n_ppn=12, walltime='6:00:00')
  assert all((tabfile, outdir))
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

  # Only re-create varlist and numpy matrix if they do not yet exist in same directory.
  varlist_fname = os.path.join(outdir, os.path.basename(tabfile) + ".varlist.txt")
  npy_fname = os.path.join(outdir, "%s.npy" % tabfile)
  if os.path.exists(varlist_fname) and os.path.exists(npy_fname):
    print "Both %s and %s exist, do not recreate varlist and numpy masked matrix files." % \
        (varlist_fname, npy_fname)
    # load numpy matrix to get its size
    n = np.size(ma.load(npy_fname), 0)
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

  # Move npy_fname to workdir
  work_npy_fname = os.path.join(work_dir, os.path.basename(npy_fname))
  if not os.path.exists(work_npy_fname):
    shutil.move(npy_fname, work_npy_fname)
    print "Moved %s to %s." % (npy_fname, work_npy_fname)

  # dispatch jobs in a loop
  num_pairs = int(n * (n-1) / 2) # no diagonal: n choose 2

  # Create new dispatch script file.o
  dispatch_script_fname = make_script_name(work_dir, work_npy_fname, "dispatch_pairwise")
  print "Creating batch script '%s'..." % dispatch_script_fname
  fp = open(dispatch_script_fname, 'w')

  # Write jobs to dispatch script in a list.
  offset = start_offset
  t = tstamp()
  while offset < num_pairs:
    # BATCH_CMD = "python %(script_path)s/batch_pairwise.py npyfile=%(npyfile)s offset=%(offset)d k=%(k)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"
    for function in all_functions:
      cmd = BATCH_CMD % {
        'script_path': os.path.dirname(os.path.realpath(__file__)),
        'npyfile': work_npy_fname,
        'offset': offset,
        'k': k,
        'n': n,
        'work_dir': outdirs[function],
        'stdout_fname': os.path.join(work_dir, "log_%s_%s_%s.out" % (jobname, t, function)),
        'stderr_fname': os.path.join(work_dir, "log_%s_%s_%s.err" % (jobname, t, function)),
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
