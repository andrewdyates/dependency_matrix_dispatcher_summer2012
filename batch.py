#!/usr/bin/python
import os
import numpy as np
from scipy.stats import mstats

def make_Ms(size, names):
  d = {}
  for n in names:
    d[n] = np.zeros(size)
  return d

class Batch(object):
  MNAMES = []
  
  def __init__(self, size=1):
    """Initialize results matrices."""
    assert size is not None and size > 0
    self.size = size
    self.Matrices = make_Ms(size, self.MNAMES)
    
  def compute(self, x, y, i):
    """Compute measure of dependencies and save to results matrices."""
    assert np.size(x) == np.size(y) and i >= 0
    raise Exception, "Not Implemented."
  
  def compute_one(self, x, y):
    """Compute one value in the first matrix entry, return it."""
    self.compute(x, y, 0)
    return self.get(0)
  
  def get(self, i):
    """Return list of values computed for index i"""
    assert i >= 0
    return dict([(n, self.Matrices[n][i]) for n in self.MNAMES])
  
  def nans(self):
    """Return total number of not-a-numbers in all matrices."""
    return sum([np.sum(np.isnan(M)) for M in self.Matrices.values()])
  
  def save(self, work_dir, batchname):
    """Save result matrices to file."""
    out_names = {}
    for name, M in self.Matrices.items():
      output_fname = os.path.join(work_dir, "%s.%s.npy" % (batchname, name))
      np.save(output_fname, M)
      out_names[name] = output_fname
    return out_names

  
class PCCBatch(Batch):
  MNAMES = ["PEARSON", "PEARSON_PV"]
  HAS_NEG = True
  def compute(self, x, y, i):
    assert np.size(x) == np.size(y) and i >= 0
    self.Matrices["PEARSON"][i], self.Matrices["PEARSON_PV"][i] = mstats.pearsonr(x,y)

class CovBatch(Batch):
  MNAMES = ["COVARIANCE", "STD_PRODUCT", "MIN_STD"]
  HAS_NEG = True
  def compute(self, x, y, i):
    assert np.size(x) == np.size(y) and i >= 0
    self.Matrices["COVARIANCE"][i] = np.cov(x,y)[0,1]
    self.Matrices["MIN_STD"][i] = min(np.std(x),np.std(y))
    self.Matrices["STD_PRODUCT"][i] = np.std(x) * np.std(y)
    
class SpearmanBatch(Batch):
  MNAMES = ["SPEARMAN", "SPEARMAN_PV"]
  HAS_NEG = True
  def compute(self, x, y, i):
    assert np.size(x) == np.size(y) and i >= 0
    self.Matrices["SPEARMAN"][i], self.Matrices["SPEARMAN_PV"][i] = mstats.spearmanr(x,y)

class EuclideanBatch(Batch):
  MNAMES = ["EUCLIDEAN"]
  def compute(self,x,y,i):
    assert np.size(x) == np.size(y) and i >= 0
    q=x-y
    self.Matrices["EUCLIDEAN"][i] = np.sqrt((q*q.T).sum())

#  'kendalltau': lambda x,y: mstats.kendalltau(x,y)[0],
class KendallBatch(Batch):
  MNAMES = ["KENDALL", "KENDALL_PV"]
  def compute(self,x,y,i):
    assert np.size(x) == np.size(y) and i >= 0
    k, p = mstats.kendalltau(x,y)
    self.Matrices["KENDALL"][i] = k
    self.Matrices["KENDALL_PV"][i] = p

