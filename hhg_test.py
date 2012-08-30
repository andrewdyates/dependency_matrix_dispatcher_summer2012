from hhg_batch import *
import numpy as np
F = HHGBatch()
for i in xrange(10, 1000):
  x = np.arange(i)
  print F.compute_one(x, x)
