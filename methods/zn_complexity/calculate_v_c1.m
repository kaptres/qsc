% author: mferhata
function calculate_v_c1 (nam, num)
    load ([ '~/SHREC_collection1/' char(nam) sprintf('/%03d',num) '.mat']);
    S   = trim_volume (S);
    dt  = bwdist (~S, 'ch');
    v   = spe_Linfty_iterative_mex_3D (S, max(dt(:)), .49, 1e-7, 200, analytical_solution(S,1,dt));
    save (['~/DATA/SHREC_collection1/' char(nam) '_Mar26/' sprintf('%03d',num) '_v.mat'], 'v');
end
