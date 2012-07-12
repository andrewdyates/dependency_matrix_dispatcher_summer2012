#!/usr/bin/python
"""Dispatch one-per-variable MINE.jar jobs using parallel-command-processor.

EXAMPLE USE:
  python dispatch_minejar.py tabfile=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab
"""
import os, sys
import subprocess
import shutil
import datetime

CMD = "python %(scriptfile) tabfile=%(tabfile)s offset=%(offset)d"
WORK_DIR = "/fs/lustre/osu6683"
MINEJAR_FILE = "/fs/lustre/osu6683/MINE.jar"
#TMP_DIR = os.environ['TMPDIR']

BATCH_CMD = "python %(script_path)s/batch_minejar.py tabfile=%(tabfile)s offset=%(offset)d minejar_file=%(minejar_file)s work_dir=%(work_dir) >> %(stdout_fname)s 2>> %(stderr_fname)s"

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

def dispatch(tabfile=None, n_nodes=2, n_ppn=12, walltime='6:00:00', \
         work_dir=WORK_DIR, minejar_file=MINEJAR_FILE, dry=False, start_offset=0, \
         jobname=None):
  """Create parallel-command-processor script and start paralell job.
  """
  assert tabfile
  n_nodes, n_ppn, start_offset = map(int, (n_nodes, n_ppn, start_offset))
  if type(dry) == str and dry.lower() in ('false', 'f', 'none', ''):
    dry=False
  if jobname is None: 
    jobname = os.path.basename(tabfile)

  # Copy tabfile to WORK_DIR.
  tab_basename = os.path.basename(tabfile)
  work_tabfile = os.path.join(work_dir, tab_basename)
  shutil.copyfile(tabfile, work_tabfile)
  print "Copied %s to %s." % (tabfile, work_tabfile)
  
  # Count number of variables in .tab file.
  n = len([f for f in open(work_tabfile) if f[0] not in ('\n', '#')])
  print "Counted %d variables in tabfile %s." % (n, work_tabfile)
  assert start_offset < n

  # Create new dispatch script file.
  script_path = os.path.dirname(os.path.realpath(__file__)) # in same path as this script
  tstamp = datetime.datetime.isoformat(datetime.datetime.now())
  tmp_script_name = "tmp_script_%s_%s_%s.sh" % (tab_basename, "dispatch_minejar.py", tstamp)
  dispatch_script_fname = os.path.join(work_dir, tmp_script_name)
  print "Creating batch script '%s'..." % dispatch_script_fname
  fp = open(dispatch_script_fname, 'w')

  # Write jobs to dispatch script in a list.
  offset = start_offset
  while offset < n:
  #    BATCH_CMD = "python %(script_path)s/batch_minejar.py tabfile=%(tabfile)s offset=%(offset)d minejar_file=%(minejar_file)s work_dir=%(work_dir) >> %(stdout_fname)s 2>> %(stderr_fname)s"
    cmd = CMD % {
      'script_path': script_path,
      'tabfile': work_tabfile,
      'offset': offset, 
      'minejar_file': minejar_file,
      'work_dir': work_dir,
      'stdout_fname': os.path.join(work_dir, "log_"+jobname+"_"+tstamp+".out"),
      'stderr_fname': os.path.join(work_dir, "log_"+jobname+"_"+tstamp+".err"),
      }
    fp.write(cmd); fp.write('\n')
    offset += 1
  fp.close()

  # Submit job script.
  qsub_script = QSUB_TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, 'dispatch_script': dispatch_script_fname}

  print qsub_script
  if not dry:
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE)
    p.communicate(input=qsub_script)
    p.stdin.close()
    print "Batch job submitted."
  else:
    print "Dry run."

  
if __name__ == "__main__":
  dispatch(**dict([s.split('=') for s in sys.argv[1:]]))
