#!/usr/bin/python
"""Dispatch one-per-variable MINE.jar jobs using parallel-command-processor.

EXAMPLE USE:
  python dispatch_minejar.py tabfile=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab
"""
import os, sys
import shutil
from util import *


MINEJAR_FILE = "/fs/lustre/osu6683/MINE.jar"
#TMP_DIR = os.environ['TMPDIR']

BATCH_CMD = "time python %(script_path)s/batch_minejar.py tabfile=%(tabfile)s offset=%(offset)d minejar_file=%(minejar_file)s work_dir=%(work_dir)s >> %(stdout_fname)s 2>> %(stderr_fname)s"

def dispatch_minejar(tabfile=None, n_nodes=2, n_ppn=12, walltime='6:00:00', \
         work_dir=WORK_DIR, minejar_file=MINEJAR_FILE, dry=False, start_offset=0, \
         jobname=None):
  """Create parallel-command-processor script and start paralell job."""
  assert tabfile
  n_nodes, n_ppn, start_offset = map(int, (n_nodes, n_ppn, start_offset))
  dry = is_dry(dry)

  if jobname is None: 
    jobname = os.path.basename(tabfile)

  # Copy tabfile to WORK_DIR.
  tab_basename = os.path.basename(tabfile)
  work_tabfile = os.path.join(work_dir, tab_basename)
  shutil.copyfile(tabfile, work_tabfile)
  print "Copied %s to %s." % (tabfile, work_tabfile)
  
  # Count number of variables in .tab file.
  n = count_tab_rows(work_tabfile)
  print "Counted %d variables in tabfile %s." % (n, work_tabfile)
  assert start_offset < n

  # Create new dispatch script file.
  dispatch_script_fname = make_script_name(work_dir, tab_basename, "dispatch_minejar.py")
  print "Creating batch script '%s'..." % dispatch_script_fname
  fp = open(dispatch_script_fname, 'w')

  # Write jobs to dispatch script in a list.
  offset = start_offset
  t = tstamp()
  while offset < n:
  #    BATCH_CMD = "python %(script_path)s/batch_minejar.py tabfile=%(tabfile)s offset=%(offset)d minejar_file=%(minejar_file)s work_dir=%(work_dir) >> %(stdout_fname)s 2>> %(stderr_fname)s"
    cmd = BATCH_CMD % {
      'script_path': os.path.dirname(os.path.realpath(__file__)) ,
      'tabfile': work_tabfile,
      'offset': offset, 
      'minejar_file': minejar_file,
      'work_dir': work_dir,
      'stdout_fname': os.path.join(work_dir, "log_"+jobname+"_"+t+".out"),
      'stderr_fname': os.path.join(work_dir, "log_"+jobname+"_"+t+".err"),
      }
    fp.write(cmd); fp.write('\n')
    offset += 1
  fp.close()

  # Submit job script.
  qsub_script = QSUB_TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, 'dispatch_script': dispatch_script_fname, 'num_jobs': n, 'k': 0, 'num_pairs': n*(n-1)/2, 'function': 'MINE.jar'}
  print qsub_script
  if not dry:
    fork_qsub_submit(qsub_script)
    print "Batch job submitted."
  else:
    print "Dry run."

  
if __name__ == "__main__":
  dispatch_minejar(**dict([s.split('=') for s in sys.argv[1:]]))
