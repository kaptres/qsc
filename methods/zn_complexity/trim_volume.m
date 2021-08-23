% author: mferhata
function S  = trim_volume (S)
    inds    = find (S==1);
    [x,y,z] = ind2sub (size(S), inds);
    S       = S(min(x):max(x), min(y):max(y), min(z):max(z));
    S       = padarray (S, [1 1 1]);
end
