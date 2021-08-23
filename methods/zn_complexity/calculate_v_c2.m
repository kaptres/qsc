% author: mferhata
function calculate_v_c2 (nam, num, init)
    OBJ = read_wobj (['~/collection2/' char(nam) sprintf('/%02d',num) '.obj']);
    FV  = OBJ2FV (OBJ);
    S   = voxelize_FV (FV);
    S   = trim_volume (S);
    dt  = bwdist (~S, 'ch');
    if exist ('init', 'var')
        v   = spe_Linfty_iterative_mex_3D (S, max(dt(:)), .49, 1e-7, 100, init);
    else
        v   = spe_Linfty_iterative_mex_3D (S, max(dt(:)), .49, 1e-7, 600, analytical_solution(S,1,dt));
    end
    save (['~/DATA/SHREC_collection2/' char(nam) '_Mar26/' sprintf('%02d',num) '_v.mat'], 'v');
end
