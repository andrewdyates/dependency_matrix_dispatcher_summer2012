import numpy as np
import dcor as py
import dcor_cpy as cy
A = np.random.rand(1000)
B = np.random.rand(1000)
print py.dcov_all(A,B)
print cy.dcov_all(A,B)
