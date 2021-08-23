% Author: Murat Genctav
% The corresponding paper:
% Genctav, M., Genctav, A. & Tari, S. 
% NonLocal via Local-NonLinear via Linear: A New Part-coding Distance Field via Screened Poisson Equation. 
% J Math Imaging Vis 55, 242-252 (2016). https://doi.org/10.1007/s10851-015-0614-8
function stepDist = parcellinStepDistance(volume,STEP,NSCALES)

[WS,~,~,~] = parcellin3Dfield(volume,STEP,NSCALES);

i = 0;
emptyPositiveSet = false;
w = zeros(size(volume));
while ~emptyPositiveSet && i < NSCALES
	w(volume == 1) = WS(:,i+1);
    if i==0
		stepDist = uint32(w<0);
	else
		stepDist = stepDist + uint32(w<0);
    end
    if ~any(any(any(w>0)))
        emptyPositiveSet = true;
    end
	i = i+1;
end

end