The files in this folder reproduce the results provided as Heatmaps.
The generated Heatmaps are saved in the calling directory as PNG files with hard-coded reference information.

'mcorr_collection1.m' creates 4 files under the calling directory:
    cp_corr.png     -> Fig. 7 top - the leftmost image
    sp_corr.png     -> Fig. 7 top - the 2nd-from-left image
    cm_corr.png     -> Fig. 7 top - the 2nd-from-right image
    sm_corr.png     -> Fig. 7 top - the rightmost image

'mcorr_collection2.m' creates 2 files under the calling directory:
    c21_corr.png    -> Fig. 7 bottom - the leftmost image
    c22_corr.png    -> Fig. 7 bottom - the 2nd-from-left image

'mcorr_collection3.m' creates 2 files under the calling directory:
    c3_cat_corr.png -> Fig. 7 bottom - the 2nd-from-right image
    c3_all_corr.png -> Fig. 7 bottom - the rightmost image

------------------------------------------------------------------------------------------

sample_run.m:                       the one-call script to generate the Heatmaps.

parse_submission2.m:                parses the results-file indicated by the variable 'filename'.
load_reverseds.m:                   reverses the signs of the methods that measure 'simplicity'.
pairwise_correlation_collectionX.m: computes Kendall's tau between two results-files on Collection X.
mcorr_collectionX.m:                computes Kendall's tau between all results-files pairs on Collection X.
                                    plots and saves the Heatmaps.

------------------------------------------------------------------------------------------
