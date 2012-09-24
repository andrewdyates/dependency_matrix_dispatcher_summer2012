dependency_matrix_dispatcher
============================
Compute all-pair dependency matrices in parallel.
All matrices are assumed to be row vectors, column dimensions.

TODO:
* move dependency computations into own module
* merge dispatches into single function
* include compilation in dispatch


*See individual file docstrings for use details.*

####Computes####
* Pearson's Correlation (R or PCC)
* Spearman's Correlation (rho)
* Distance Correlation (dCOR)
* MINE (using the minepy module at http://minepy.sourceforge.net/)
* HHG (Heller Heller Gorfine 2012 statistic)

Workflow
--------

1. `dispatch_*.py` reads a data matrix, divides work into many executions of `batch_*.py`,
and dispatches these jobs using `qsub` to compute this work in parallel.
2. `batch_*.py` computes a set of vector pair dependencies and saves the result as an
enumerated binary numpy (*.npy) matrix.
3. `compile_*.py` reads all the *.npy matrices saved from each completed `batch_*.py`
execution and compiles them all into a single .npy matrix


###Compute n choose 2 pairs of dependencies.###

1. dispatch_pairwise.py
2. batch_pairwise.py
3. compile_pairwise.py

###Compute n*m pairs of dependencies between two matrices.###

1. dispatch_dual_pairwise.py
2. batch_dual_pairwise.py
3. compile_dual_pairwise.py

