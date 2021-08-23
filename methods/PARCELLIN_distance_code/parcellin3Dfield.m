% Author: Murat Genctav
% Corresponding paper:
% Genctav, M., Genctav, A. & Tari, S. 
% NonLocal via Local-NonLinear via Linear: A New Part-coding Distance Field via Screened Poisson Equation. 
% J Math Imaging Vis 55, 242-252 (2016). https://doi.org/10.1007/s10851-015-0614-8
function [WS,timeSolve,timeSetup,timeDT] = parcellin3Dfield(shapeMask,scaleRes,numScales)
% PARCELLIN3DFIELD computes multiple PARCELLIN fields inside the shape for
% varying s parameters from 0 to "numScales*scaleRes" with step size "scaleRes".
% INPUT:
%   "shapeMask" is a 3D matrix that stores shape indicator function.
%   Attains 1 inside, and 0 outside the shape.
%   "scaleRes" determines the sampling resolution of the scales at which the
%   field is to be computed.
%	"numScales" determines the number of different scales.
% COMPUTATION:
%   L(ws)-alpha*ws= -dist
%   where
%   L denotes the Laplacian,
%   alpha may be selected as 1/(number of nodes)^2 inside the shape domain,
%   dist may be the Euclidean distance function calculated inside the shape.

% ensure the shape domain does not touch the frame boundary
shapeMask = padarray(shapeMask,[1 1 1]);
% compute screening parameter
alpha = 1/nnz(shapeMask)^2;

% compute EDT
%fprintf('\t\tComputing normalized Euclidean Distance Transform...');
tic
EDT = bwdist(~shapeMask);
EDT = double(EDT/max(EDT(:)));
timeDT = toc;
%fprintf('Done in %f seconds\n',t);

% set up screened Poisson system matrix with multiple RHS
%fprintf('\t\tSetting up screened Poisson system matrix and multiple RHSs...');
tic
[A,B] = buildScreenedPoissonSystemMex(EDT,alpha,scaleRes,numScales);
timeSetup = toc;
%fprintf('Done in %f seconds\n',t);

% solve the equation using direct Cholesky solver
%fprintf('\t\tRunning direct Cholesky solver with multiple RHSs...\n');
tic
WS = -A\B;
timeSolve = toc;
%fprintf('\t\tDone in %f seconds\n',st);
end