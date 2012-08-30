from batch import Batch
# R and Rpy2 must be installed.
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import numpy as np

class HHGBatch(Batch):
  MNAMES = ["SUM_CHI", "SUM_LR", "MAX_CHI", "MAX_LR"]

  def __init__(self, size):
    # Load HHG library from R installation.
    r('library("HHG2x2")')
    super(HHGBatch, self).__init__(size)

  def compute(self, x, y, i):
    assert type(int(i)) == int
    assert np.size(x) == np.size(y) and i >= 0
    n = np.size(x)
    robjects.globalenv["x"] = x
    robjects.globalenv["y"] = y
    r('Dx = as.matrix(dist((x),diag=TRUE,upper=TRUE))')
    r('Dy = as.matrix(dist((y),diag=TRUE,upper=TRUE))')
    HHG = r('myHHG(Dx,Dy)')
    v = n*(n-2)*(n-3) # max sum_chi value
    lg = np.log(2)
    self.Matrices["SUM_CHI"][i] = float(HHG.rx('sum_chisquared')[0][0]) / v
    self.Matrices["SUM_LR"][i]  = float(HHG.rx('sum_lr')[0][0]) / v / lg
    self.Matrices["MAX_CHI"][i] = float(HHG.rx('max_chisquared')[0][0]) / (n-2)
    self.Matrices["MAX_LR"][i]  = float(HHG.rx('max_lr')[0][0]) / (n-2) / lg
