dependency_matrix_dispatcher
============================
Compute all-pair dependency matrices in parallel.
All matrices are assumed to be row vectors, column dimensions.

*See individual file docstrings for use details.*

####Computes####
* Pearson's Correlation (R or PCC)
* Spearman's Correlation (rho)
* Distance Correlation (dCOR)
* MINE (using the minepy module at http://minepy.sourceforge.net/)

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

Alternate batch_*.py for dispatch_pairwise.py
* batch_minepy_pairwise.py: Save 4 MINE statistics (MIC, MEV, MCN, MAS) computed using cmine Python wrapper

###Compute n choose 2 pairs of dependencies using Reshef et al. MINE.jar###
_depreciated_

1. dispatch_minejar.py
2. batch_minejar.py
3. compile_minejar.py

###Compute n*m pairs of dependencies between two matrices.###

1. dispatch_dual_pairwise.py
2. batch_dual_pairwise.py
3. compile_dual_pairwise.py



