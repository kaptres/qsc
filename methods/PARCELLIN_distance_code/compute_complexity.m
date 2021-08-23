function val = compute_complexity(volume)

STEP = 0.01;
NSCALES = 70;
stepDist = parcellinStepDistance(volume,STEP,NSCALES);

se = strel('cube',3);
volBoundary = volume-imerode(volume,se);

% get step dist values of boundary elements
boundDistValues = stepDist(volBoundary == 1);

% compute histogram
binedges = (0.5:1:NSCALES+0.5);
histCount = histcounts(boundDistValues,binedges);

% normalize histogram counts
p = histCount / length(boundDistValues);

% compute entropy as the complexity
pp = p(p>0);
val = -sum(pp.*log2(pp));

end