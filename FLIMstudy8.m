%% Under sample and SNR
% Load
load("sampledata_dlight.mat");



%% Simulation parameters
% 21 st data point as T = 0
t0 = t(21);

% Photons to sim
photons2sim = [25 50 100 200 500 1000 5000 10000 100000 350000];
nsims = length(photons2sim);

% Iterations
niter = 100;

% Scale deconv
deconv2 = round(deconv*19000);

% Ground truth
fprintf('Ground truth Tm: %0.3f ns\n', sum(raw .* t) / sum(raw) - t0);
fprintf('Ground truth IEM: %0.3f ns\n', sum(deconv2) / max(deconv2) * t(2));

%% Simulation Tm
% t matrix
tmat = t * ones(1,niter);

% Initialize
simdata_tm = zeros(niter, nsims);

hwait = waitbar(0);
for i = 1 : nsims
    tic;
    waitbar(i/nsims, hwait, sprintf('%i/%i', i, nsims))
    hists_sim_tm = histresample(raw, photons2sim(i), niter);
    simdata_tm(:,i) = (sum(hists_sim_tm .* tmat, 1) / photons2sim(i))' - t0;
    fprintf('%i: %0.1f s\n', i, toc)
end
close(hwait)

%% Descriptive stats
upper_tm = prctile(simdata_tm, 97.5, 1);
lower_tm = prctile(simdata_tm, 2.5, 1);
meanval_tm = mean(simdata_tm, 1);
stddev_tm = std(simdata_tm);

%% Simulation IEM
% t matrix
tmat = t * ones(1,niter);

% Initialize
simdata_iem = zeros(niter, nsims);

hwait = waitbar(0);
for i = 1 : nsims
    tic;
    waitbar(i/nsims, hwait, sprintf('%i/%i', i, nsims))
    hists_sim_iem = histresample(deconv2, photons2sim(i), niter);
    simdata_iem(:,i) = sum(hists_sim_iem,1) ./ max(hists_sim_iem, [], 1) * t(2);
    fprintf('%i: %0.1f s\n', i, toc)
end
close(hwait)

%% Descriptive stats
upper_iem = prctile(simdata_iem, 97.5, 1);
lower_iem = prctile(simdata_iem, 2.5, 1);
meanval_iem = mean(simdata_iem, 1);
stddev_iem = std(simdata_iem) ./ meanval_iem;