function newhist = histresample(inputhist, n, niter)
% histresample resamples histogram data with N points
% newhist = histresample(inputhist, n, niter)

if nargin < 3
    niter = 1;
end

%% Constants
% total number of datapoints
total = sum(inputhist);

% N bins
nbins = length(inputhist);

% Cumulative histogram
cumsum_his = cumsum(inputhist);

%% Sample
% New histogram
newcumhist = zeros(nbins, niter);

% random sample
inds = randi(total, [n, niter]);

% Generate new cumulative histogram
for ibin = 1 : nbins
    for iter = 1 : niter
        newcumhist(ibin, iter) = sum(inds(:,iter) <= cumsum_his(ibin));
    end
end

% New histograms
newhist = cat(1, newcumhist(1,:), diff(newcumhist, 1, 1));

end