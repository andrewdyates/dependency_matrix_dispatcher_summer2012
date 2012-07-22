#!/usr/bin/python
"""Dispatch pairwise dependency computations between two matrices.

USE EXAMPLE:

python $HOME/dependency_matrix_dispatcher/dispatch_dual_pairwise.py outdir=/fs/lustre/osu6683/gse15745 tabfile_1=$HOME/gse15745/GSE15745.GPL6104.mRNA.normed.tab tabfile_2=$HOME/gse15745/GSE15745.GPL8178.miRNA.normed.tab tabfile_1_coltitles=$HOME/gse15745/GSE15745.GSE6104.subject_titles.txt tabfile_2_coltitles=$HOME/gse15745/GSE15745.GSE8178.subject_titles.txt dry=True
"""
from __future__ import division
from util import *
import sys
import shutil
import numpy.ma as ma


BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s > %(stdout_fname)s 2> %(stderr_fname)s"

def dispatch_dual_pairwise(tabfile_1=None, tabfile_2=None, tabfile_1_coltitles=None, tabfile_2_coltitles=None, outdir=WORK_DIR, function=None, k=200000, dry=False, start_offset=0, work_dir=WORK_DIR, jobname=None, n_nodes=6, n_ppn=12, walltime='12:00:00'):
  assert tabfile_1 and tabfile_2 and tabfile_1_coltitles and tabfile_2_coltitles
  n_nodes, n_ppn, start_offset, k = map(int, (n_nodes, n_ppn, start_offset, k))
  assert k > 1 and start_offset >= 0 and n_nodes > 0 and n_ppn > 0
  if jobname is None: 
    jobname_base = "%s_vs_%s" % (os.path.basename(tabfile_1), os.path.basename(tabfile_2))
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

  npy_fname_1, n1, M1 = npy_varlist_from_tabfile(tabfile_1, outdir, overwrite=True)
  npy_fname_2, n2, M2 = npy_varlist_from_tabfile(tabfile_2, outdir, overwrite=True)

  # let the smaller matrix be matrix 1
  if n1 > n2:
    print "Matrix 1 [%s] has more rows (%d) than Matrix 2 [%s] (%d)" % \
        (npy_fname_1, n1, npy_fname_2, n2)
    npy_fname_1, npy_fname_2 = npy_fname_2, npy_fname_1
    tabfile_1_coltitles, tabfile_2_coltitles = tabfile_2_coltitles, tabfile_1_coltitles
    n1, n2 = n2, n1
    M1, M2 = M2, M1
    print "Swapped matrices 1 and 2 so that matrix 1 is smaller."

  # Print analytics
  print "Created Matrix 1: %s" % npy_fname_1
  print_matrix_stats(M1)
  print "Created Matrix 2: %s" % npy_fname_2
  print_matrix_stats(M2)
  
  # Column-Align two matrices. 
  npy_basename_1 = os.path.basename(npy_fname_1)
  npy_basename_2 = os.path.basename(npy_fname_2)
  npy_aligned_fname_1 = os.path.join(work_dir, "%s_aligned_with_%s.npy" % (npy_basename_1, npy_basename_2))
  npy_aligned_fname_2 = os.path.join(work_dir, "%s_aligned_with_%s.npy" % (npy_basename_2, npy_basename_1))
  titles1, set1, idx1 = read_samples(tabfile_1_coltitles)
  titles2, set2, idx2 = read_samples(tabfile_2_coltitles)
  cols = []
  for s in titles1:
    if s in set2:
      cols.append((idx1[s], idx2[s]))
  # Align columns in copies of matrices.
  print "Aligned %d cols from %s and %d cols from %s to %d shared columns by sample ID." % \
    (len(titles1), tabfile_1_coltitles, len(titles2), tabfile_2_coltitles, len(cols))
  m1_idxs, m2_idxs = np.array(zip(*[(x, y) for x, y in cols]))
  M1_aligned = M1[:,m1_idxs]
  M2_aligned = M2[:,m2_idxs]
  print "Created (%d x %d) M1 and (%d x %d) M2" % (n1, len(cols), n2, len(cols))
  print "M1_aligned: (%d x %d). M2_aligned: (%d x %d)." % \
    (np.size(M1_aligned,0),np.size(M1_aligned,1),np.size(M2_aligned,0),np.size(M2_aligned,1))
  print "M1_aligned"
  print_matrix_stats(M1_aligned)
  print "M2_aligned"
  print_matrix_stats(M2_aligned)
  
  # Save column aligned matrices.
  print "Saving %s and %s..." % (npy_aligned_fname_1, npy_aligned_fname_2)
  ma.dump(M1_aligned, npy_aligned_fname_1)
  ma.dump(M2_aligned, npy_aligned_fname_2)
  # Save column titles, verify that they are equal
  aligned_col_titles = []
  for x, y in cols:
    s1 = titles1[x]
    s2 = titles2[y]
    assert s1 == s2
    aligned_col_titles.append(s1)
  aligned_col_titles_fname = os.path.join(work_dir, "%s_vs_%s.coltitles.txt" % (npy_basename_1, npy_basename_2))
  fp = open(aligned_col_titles_fname, "w")
  fp.write(",".join(aligned_col_titles)); fp.write("\n")
  fp.close()
  print "Wrote %d aligned column titles at %s." % (len(aligned_col_titles), aligned_col_titles_fname)
  assert len(aligned_col_titles) == np.size(M1_aligned, 1) == np.size(M2_aligned, 1)
  assert np.size(M1_aligned, 0) <= np.size(M2_aligned, 0)
  
  # Write jobs to dispatch script in a list.
  t = tstamp()
  for function in all_functions:
    jobname = "%s_%s" % (jobname_base, function)
    # Create new dispatch script per function
    dispatch_script_fname = \
      make_script_name(work_dir, "%s_vs_%s" % (npy_basename_1, npy_basename_2), "dispatch_%s" % function)
    print "Creating batch script '%s'..." % dispatch_script_fname
    fp = open(dispatch_script_fname, 'w')
    offset = start_offset
    num_jobs = 0
    for offset in xrange(n1):
# BATCH_CMD = "time python %(script_path)s/batch_dual_pairwise.py npyfile_1=%(npyfile_1)s npyfile_2=%(npyfile_2)s offset=%(offset)d work_dir=%(work_dir)s function=%(function)s >> %(stdout_fname)s 2>> %(stderr_fname)s"
      cmd = BATCH_CMD % {
        'script_path': os.path.dirname(os.path.realpath(__file__)),
        'npyfile_1': npy_aligned_fname_1,
        'npyfile_2': npy_aligned_fname_2,
        'offset': offset,
        'work_dir': outdirs[function],
        'stdout_fname': os.path.join(outdirs[function], "log_%s_%s_%s_%d.out" % (jobname, t, function, offset)),
        'stderr_fname': os.path.join(outdirs[function], "log_%s_%s_%s_%d.err" % (jobname, t, function, offset)),
        'function': function,
      }
      fp.write(cmd); fp.write('\n')
      num_jobs += 1
    fp.close()
    # Submit job script.
    qsub_script = QSUB_TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, 'dispatch_script': dispatch_script_fname, 'num_jobs': num_jobs, 'function': function, 'k': n2, 'num_pairs': n1*n2}
    print qsub_script
    if not dry:
      fork_qsub_submit(qsub_script)
      print "Batch job submitted."
    else:
      print "Dry run."


if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  dispatch_dual_pairwise(**args)
