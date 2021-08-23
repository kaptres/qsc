% author: mferhata
% calculates and prints the E of the shape specified by
%       NAM: subcollection name,        NUM: shape no
%
% the collections are hard-coded: change collection2 to collection1 or collection3 accordingly
%
% v-fields should be calculate before calling this (see: calculate_v_c[1,2,3].m)
function Es = calculate_E (nam, num)
    load (['~/DATA/SHREC_collection2/' char(nam) '_Mar26/' sprintf('%02d',num) '_v.mat']);
    v       = v / max(v(:));
    [s,Es]  = TIP_score (v);
    save (['~/DATA/SHREC_collection2/' char(nam) '_Mar26/' sprintf('%02d',num) '_Es.mat'], 'Es');
end
