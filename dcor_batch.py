from batch import Batch
import numpy as np
# Comment cpy dcor, uncomment native python dcor for pychecker.
from dcor_cpy import *
#from dcor import *


class DcorBatch(Batch):
  MNAMES = ["DCOR", "DCOV"]

  def compute(self, x, y, i):
    assert type(int(i)) == int
    assert np.size(x) == np.size(y) and i >= 0
    dc, dr, dvx, dvy = dcov_all(x,y)
    self.Matrices["DCOR"][i] = dr
    self.Matrices["DCOV"][i] = dc
