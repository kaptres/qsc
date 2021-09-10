The files in this folder reproduce the results reported as Kendall-tau scores of a method.
The scores are to indicate the methods' performances with respect to the ground truth (GT).
This will generate a source code for a LaTeX table with hard-coded reference information and table alignment.

To acquire results call 'sample_run.m' which by default generates the results for the additively perturbed sets of Collection 1.
The other results can be acquired by uncommenting the related lines.

------------------------------------------------------------------------------------------

sample_run.m:           the one-call script to generate the tables.

parse_submission2.m:    parses the results-file indicated by the variable 'filename'.
kendall_tau.m:          manual implementation of the Kendall's tau computation.
segmentations.m:        initializes the GT for Collection 3.
collectionX_scores.m:   generates the Kendall's tau-based scores of a given method for Collection X.
all_scores.m:           sets up the related variables for the scripts 'create_colX_table.m' to use by calling the scripts 'collectionX_scores.m'.
create_colX_table.m:    creates the LaTeX tables for Collection X.

------------------------------------------------------------------------------------------
