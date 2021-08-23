% author: mferhata
function calculate_v_c3 (num)
    filepath = ['~/off/' sprintf('%d.off', num)];
    if exist (filepath) == 0
        return;
    end
    FV  = load_off_model (filepath);
    S   = voxelize_FV (FV);
    S   = trim_volume (S);
    dt  = bwdist (~S, 'ch');
    v   = spe_Linfty_iterative_mex_3D (S, max(dt(:)), .49, 1e-7, 600, analytical_solution(S,1,dt));
    save (['~/DATA/SHREC_collection3/off_Mar26/' sprintf('%03d',num) '_v.mat'], 'v');
end
