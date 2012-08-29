#!/usr/bin/python
"""As called by dispatch_minejar

OBSOLETE

python /nfs/01/osu6683/dependency_matrix_dispatcher/batch_minejar.py tabfile=/fs/lustre/osu6683/GSE15745.GPL6104.mRNA.normed.tab offset=22180 minejar_file=/fs/lustre/osu6683/MINE.jar work_dir=/fs/lustre/osu6683
"""
import sys, os
import subprocess
import shutil

try:
  TMPDIR = os.environ['TMPDIR']
except KeyError:
  TMPDIR = os.environ['HOME']

MINEJAR_CMD = "java -jar %(minejar_file)s %(tabfile)s 0"

def run_mine(tabfile=None, offset=None, minejar_file=None, work_dir=None):
  """ """
  assert None not in (tabfile, offset, minejar_file, work_dir)
  offset = int(offset)

  # Copy tabfile to tmpdir from offset.
  tabfile_basename = os.path.basename(tabfile)
  tab_tmpfile = os.path.join(TMPDIR, "%s_%d.tab" % \
                               (tabfile_basename.rpartition('.')[0], offset))
  fp_in = open(tabfile, 'r')
  fp_out = open(tab_tmpfile, 'w')
  i = 0
  n_lines = 0
  while i < offset:
    line = fp_in.next()
    if line[0] not in ('\n', '#'):
      i += 1
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
  for filename in os.listdir(TMPDIR):
    src = os.path.join(TMPDIR,filename)
    if "Results" in filename:
      dst = os.path.join(work_dir,filename)
      print "Moving %s to %s." % (src, dst)
      shutil.move(os.path.join(TMPDIR,src), dst)
    else:
      print "Deleting %s." % src
      try:
        os.remove(src)
      except OSError, e:
        print "Error:", e



if __name__ == "__main__":
  run_mine(**dict([s.split('=') for s in sys.argv[1:]]))
