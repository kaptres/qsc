% dirPath : path of the directory containing 2D views of input 3D shape
% shapeName : name of input 3D shape 

% val : complexity value of the input shape

function val = compute_complexity(dirPath, shapeName)

files = dir(fullfile(dirPath,[shapeName '*']));

entropy_all = [];

for i = 1:length(files)
    mask = imread(fullfile(dirPath, files(i).name));
    if length(size(mask)) == 3
        mask = rgb2gray(mask);
    end
    mask = mask == 0;
    [vbes,~,~,~,~] = bessel2Dfield(mask);
    edges = -1:0.01:1;
    N = histcounts(vbes(mask),edges);
    N = N/sum(N);
    entropy = -sum(N(N>0).*log2(N(N>0)));
    entropy_all = [entropy_all entropy];
end

val = mean(entropy_all);

end