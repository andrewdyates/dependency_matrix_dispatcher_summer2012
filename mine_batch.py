#!/usr/bin/python
# minepy must be installed. See http://mpba.fbk.eu/cmine
import minepy

from batch import *

class MINEBatch(Batch):
  MNAMES = ["MIC", "MAS", "MEV", "MCN"]

  def __init__(self, alpha=0.6, c=15, *args, **kwds):
    self.mine = minepy.MINE(alpha=alpha, c=c)
    super(MINEBatch, self).__init__(*args, **kwds)
    
  def compute(self, x, y, i):
    assert type(int(i)) == int
    assert np.size(x) == np.size(y) and i >= 0
    self.mine.score(x, y)
    self.Matrices['MIC'][i] = self.mine.mic()
    self.Matrices['MAS'][i] = self.mine.mas()
    self.Matrices['MEV'][i] = self.mine.mev()
    self.Matrices['MCN'][i] = self.mine.mev()

