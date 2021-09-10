The files in this folder reproduce the results provided as Heatmaps.
The generated Heatmaps are saved in the calling directory as PNG files with hard-coded reference information.

------------------------------------------------------------------------------------------

sample_run.m:                       the one-call script to generate the Heatmaps.

parse_submission2.m:                parses the results-file indicated by the variable 'filename'.
load_reverseds.m:                   reverses the signs of the methods that measure 'simplicity'.
pairwise_correlation_collectionX.m: computes Kendall's tau between two results-files on Collection X.
mcorr_collectionX.m:                computes Kendall's tau between all results-files pairs on Collection X.
                                    plots and saves the Heatmaps.

------------------------------------------------------------------------------------------
