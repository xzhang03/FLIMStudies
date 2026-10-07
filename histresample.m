function newhist = histresample(inputhist, n, niter)
% histresample resamples histogram data with N points
% newhist = histresample(inputhist, n, niter)

if nargin < 3
    niter = 1;
end

%% Fast path: compiled multinomial sampler (mex -O histresample_mex.cpp)
% Same distribution, cost independent of n. Seeded from MATLAB's rng.
if exist('histresample_mex', 'file') == 3
    newhist = histresample_mex(double(inputhist), n, niter, randi(2^32-1));
    return
end

%% Constants
% total number of datapoints
total = sum(inputhist);

% N bins
nbins = length(inputhist);

% Lookup table: datapoint index -> bin it belongs to
bin_lut = repelem((1:nbins)', inputhist(:));

%% Sample
% New histogram
newhist = zeros(nbins, niter);

for iter = 1 : niter
    % random sample (with replacement)
    inds = randi(total, [n, 1]);

    % Count samples per bin
    newhist(:, iter) = accumarray(bin_lut(inds), 1, [nbins, 1]);
end

end