% Asli Genctav, Sibel Tari, Discrepancy: Local/Global Shape Characterization with a Roundness Bias. 
% J Math Imaging Vis 61, 160-171 (2019). https://doi.org/10.1007/s10851-018-0851-8
function [vbes,v,bes,R,maxR] = bessel2Dfield(shapeMask)
% ensure shape domain does not touch the frame boundary
shapeMask = padarray(logical(shapeMask),[1 1]);

edt = double(bwdist(~shapeMask));
maxR = max(max(edt));

% v field
rho = maxR;
v = ones(size(shapeMask));
v(shapeMask) = vfield2D_cholesky(shapeMask,rho);

R = zeros(size(shapeMask));
R(shapeMask) = maxR - edt(shapeMask);
alpha = 1 / rho;
bes = ones(size(shapeMask));
bes(shapeMask) = besseli(0,alpha*R(shapeMask)) / besseli(0,alpha*maxR);

vbes = v - bes;

vbes = vbes(2:end-1,2:end-1);
v = v(2:end-1,2:end-1);
bes = bes(2:end-1,2:end-1);
R = R(2:end-1,2:end-1);
end