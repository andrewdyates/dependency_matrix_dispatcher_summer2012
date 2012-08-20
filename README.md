dependency_matrix_dispatcher
============================
Compute all-pair dependency matrices in parallel.
All matrices are assumed to be row vectors, column dimensions.

*See individual file docstrings for use details.*

####Natively computes####
* Pearson's Correlation
* Spearman's Correlation
* Distance Correlation
* _Support for other measures of vector pair dependencies_
* _Also supports row-vs-all MINE computation using Java wrapper._

Workflow
--------

1. `dispatch_*.py` reads a data matrix, divides work into many executions of `batch_*.py`,
and dispatches these jobs using `qsub` to compute this work in parallel.
2. `batch_*.py' computes a set of vector pair dependencies and saves the result as an
enumerated binary numpy (*.npy) matrix.
3. `compile_*.py' reads all the *.npy matrices saved from each completed `batch_*.py`
execution and compiles them all into a single .npy matrix


###Compute n choose 2 pairs of dependencies.###

1. dispatch_pairwise.py
2. batch_pairwise.py
3. compile_pairwise.py

###Compute n choose 2 pairs of dependencies using Reshef et al. MINE.jar###

1. dispatch_minejar.py
2. batch_minejar.py
3. compile_minejar.py

###Compute n*m pairs of dependencies between two matrices.###

1. dispatch_dual_pairwise.py
2. batch_dual_pairwise.py
3. compile_dual_pairwise.py



