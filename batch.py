def make_Ms(size, names):
  d = {}
  for n in names:
    d[n] = np.zeros(size)
  return d

class Batch(object):
  
  def __init__(self, size):
    """Initialize results matrices."""
    self.Matrices = {}
    
  def compute(self, x, y, i):
    """Compute measure of dependencies and save to results matrices."""
    raise NotImplemented
  
  def nans(self):
    """Return total number of not-a-numbers in all matrices."""
    return sum([np.sum(np.isnan(M)) for M in self.Matrices.values()])
  
  def save(self, work_dir, batchname):
    """Save result matrices to file."""
    out_names = []
    for name, M in self.Matrices:
      output_fname = os.path.join(work_dir, "%s.%s.npy" % (batchname, name))
      np.save(output_fname, M)
      out_names.append(output_fname)
    return out_names
      

class MINEBatch(Batch):

  def __init__(self, size=None, alpha=0.6, c=15, *args, **kwds):
    assert size is not None and size > 0
    self.mine = minepy.MINE(alpha=alpha, c=c)
    self.Matrices = make_Ms(size, ["MIC", "MAS", "MEV", "MCN"])
    
  def compute(self, i, x, y):
    self.mine.score(x, y)
    self.Matrices['MIC'][i] = mine.mic()
    self.Matrices['MAS'][i] = mine.mas()
    self.Matrices['MEV'][i] = mine.mev()
    self.Matrices['MCN'][i] = mine.mev()

    
# TODO
class HHGBatch(Batch):

  def __init__(self, size=None, alpha=0.6, c=15, *args, **kwds):
    assert size is not None and size > 0
    self.mine = minepy.MINE(alpha=alpha, c=c)
    self.Matrices = make_Ms(size, ["MIC", "MAS", "MEV", "MCN"])
    
  def compute(self, i, x, y):
    self.mine.score(x, y)
    self.Matrices['MIC'][i] = mine.mic()
    self.Matrices['MAS'][i] = mine.mas()
    self.Matrices['MEV'][i] = mine.mev()
    self.Matrices['MCN'][i] = mine.mev()

