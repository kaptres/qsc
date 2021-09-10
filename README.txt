The files in this repository replicate the results presented in the paper "SHREC’21: Quantifying shape complexity".
The SHREC track's webpage: http://users.metu.edu.tr/ferhata/

The file structure is as follows:
    -> gt_kendall_tau:      The files in this folder reproduce the results reported as Kendall-tau scores of a method.
                            The scores are to indicate the methods' performances with respect to the ground truth (GT).
                            This will generate a source code for a LaTeX table with hard-coded reference information and table alignment.

    -> heatmap_generation:  The files in this folder reproduce the results provided as Heatmaps.
                            The generated Heatmaps are saved in the calling directory as PNG files with hard-coded reference information.
                            
    -> methods:             The files in this folder implement the methods that have participated in the work.
                            As the participants have submitted their end results, this folder can be considered to be there for completeness.
                            (The submitted results can be found both in 'gt_kendall_tau/results/' and 'heatmap_generation/results/'.)

To acquire sample results, provide your scores in the respective 'results' folders of gt_kendall_tau/ and heatmap_generation/, and run the respective 'sample_run.m's within MATLAB.
By default, the 'sample_run.m's produce the results reported in the paper.
