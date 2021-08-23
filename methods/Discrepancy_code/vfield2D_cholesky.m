% Author Murat Genctav muratgenctav@gmail.com
function v = vfield2D_cholesky(shapeMask,rho)
% L(v) - (1/rho^2) * v = 0
% where v = 1 at boundary

% ensure the shape domain does not touch the frame boundary
shapeMask = padarray(shapeMask,[1 1]);

% compute screening parameter
alpha = 1/rho^2;

% set up screened Poisson system matrix with multiple RHS
[A,b] = buildScreenedPoissonSystemMex(double(shapeMask),alpha);

% solve the equation using direct Cholesky solver
v = A\b;

end