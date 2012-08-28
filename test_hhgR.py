# R and Rpy2 must be installed.
# See r.matrix handling at:
# http://rpy.sourceforge.net/rpy2/doc-2.1/html/robjects.html
import numpy as np
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

# try to run R script
#data = np.random.random((10,10))
#r.heatmap(data)

# make sure that HHG2x2 has been correctly installed and compiled; see HHG_R/README.md
x = np.array([1,2,3,4,5,6,7,8,10,12,14,100,1002,1003])
y = np.array([1,2,3,4,5,6,7,1,10,12,18,100,1002,1003])
n = np.size(x)

r('library("HHG2x2", lib.loc="HHG_R")')
robjects.globalenv["x"] = x
robjects.globalenv["y"] = y

r('Dx = as.matrix(dist((x),diag=TRUE,upper=TRUE))')
r('Dy = as.matrix(dist((y),diag=TRUE,upper=TRUE))')
print r('myHHG(Dx,Dy)')

#print X[0], X.rx(1, 1)[0] 
#print X[1], X.rx(2, 1)[0]
#print X[2], X.rx(1, 2)[0]
#print "columns:", X.ncol
#print "rows:", X.nrow
# This is not elegant, but it works. Look for R Matrix row/column selection in the future.


# This works.
# r('library("HHG2x2", lib.loc="HHG_R")')
# r('X = datagenCircle(50);')
# r('Dx = as.matrix(dist((X[1,]),diag=TRUE,upper=TRUE))')
# r('Dy = as.matrix(dist((X[2,]),diag=TRUE,upper=TRUE))')
# print r('myHHG(Dx,Dy)')
