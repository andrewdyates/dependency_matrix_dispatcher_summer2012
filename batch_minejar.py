#!/usr/bin/python
"""As called by dispatch_minejar"""
import sys, os
import subprocess

BATCH_CMD = "python %(script_path)s/batch_minejar.py tabfile=%(tabfile)s offset=%(offset)d minejar_file=%(minejar_file)s work_dir=%(work_dir)s >> %(stdout_fname)s 2>> %(stderr_fname)s"
TMPDIR = os.environ['TMPDIR']

MINEJAR_CMD = "java -jar %(minejar_file)s %(tabfile)s 0"

def run_mine(script_path=None, tabfile=None, offset=None, minejar_file=None, work_dir=None, dry=False):
  """ """
  assert None not in (script_path, tabfile, offset, minejar_file, work_dir)
  offset = int(offset)
  if type(dry) == str and dry.lower() in ('false', 'f', 'none', ''):
    dry=False
  else:
    dry= True

  # Copy tabfile to tmpdir from offset.
  tab_tmpfile = os.path.join(TMPDIR, "%s_%d.tab" % (tabfile.rpartition('.')[0], offset))
  fp_in = open(tabfile, 'r')
  fp_out = open(tab_tmpfile, 'w')
  i = 0
  n_lines = 0
  while i < offset:
    fp_in.next()
  for line in fp_in:
    fp_out.write(line)
    n_lines += 1
  fp_in.close(); fp_out.close()
  print "Wrote %d lines from %s to %s." % (n_lines, tabfile, tab_tmpfile)

  # Execute java program.
  cmd = MINEJAR_CMD % {'minejar_file': minejar_file, 'tabfile': tab_tmpfile}
  print cmd
  subprocess.call(cmd, shell=True)

  # Clean up.
  for filename in TMPDIR:
    src = os.path.join(TMPDIR,filename)
    if "Results" in filename:
      dst = os.path.join(work_dir,filename)
      os.rename(src, dst)
      print "Moved %s to %s." % (src, dst)
    else:
      os.remove(src)
      print "Deleted %s." % src


if __name__ == "__main__":
  run_mine(**dict([s.split('=') for s in sys.argv[1:]]))
