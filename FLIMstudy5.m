%% Study other ways to calculate IEM with higher SNR
% The goal of this simulation is to find methods to calculate IEM with
% higher SNR. This simulation treats intensity and lifetime as two 
% independent properties.

% Simulate a sensor with two states: a (500 ps tau) and b (2000 ps tau)
% It will work out that the area if amplitudes are set to 1, the area under
% curve will also be 500 for a and 2000 for b.

% Terminology
% Amplitude: the scaler factor before the expontial term
% Intensity: the area under curve 
% Tau: time constant of the exponential term
% [Intensity] = [Amplitude] * exp(t / tau);
% Tm: first moment ("mean arrival time")
% IEM: Ratio of Intensity/Amplitude
% IEM99: Ratio of Intensity/99th percentile amplitude (instead of max)

clear
tbin = 50;
t = (1 : tbin: 10000)'; % in units of ps
atau = 500; % asum is also tau_a, in ps
btau = 2000; % sum is also tau_b, in ps
a = exp(-t/atau); % Lifetime distribution of a
b = exp(-t/btau); % Lifetime distribution of b

fprintf('Sum a-%0.1f b-%0.1f\n', sum(a), sum(b));
asum = sum(a);
bsum = sum(b);
a2 = a/asum*100;
b2 = b/bsum*100;
fprintf('Sum a2-%0.1f b2-%0.1f\n', sum(a2), sum(b2));

figure
subplot(1,2,1)
plot(t,a,t,b);
xlabel('Time (ps)')
legend({'A (tau = 500 ps)'; 'B (tau = 2000 ps)'})
title('Equal amplitude (not equal intensity)')

subplot(1,2,2)
plot(t,a/asum*100,t,b/bsum*100);
xlabel('Time (ps)')
legend({'A (tau = 500 ps)'; 'B (tau = 2000 ps)'})
title('Equal intensity (not equal amplitude)')

%% Generate data
% Experimental parameters
dff = 4;
l = 100;

% IEM99
percentile = 99;

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;
C = a2 * i1 + (b2 * i2) * dff;

% Ground truth
tm = sum(C .* (t * ones(1, l*2)), 1) ./ sum(C);
iem = sum(C) * tbin ./ max(C);
iem99 = sum(C) * tbin ./ prctile(C, percentile);

% Figure
figure
plot([tm', iem', iem99'])
legend({'Tm', 'IEM', 'IEM99'})
title('Noiseless data')

%% Test different metrics with noisy data
% Conclusion is that IEM99 is somewhat more robust to noise than IEM.
% Presumably due to the votality of extremas.

% Experimental parameters
dff = 2;
l = 100;

% IEM99
percentile = 99;

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;
C = a2 * i1 + (b2 * i2) * dff;

% Generate noise
sigma = 0.3;
N = randn(size(C)) * sigma + 1;
N(N < 0) = 0;
CN = C .* N;

% New lifetimes
tm_noise = sum(CN .* (t * ones(1, l*2)), 1) ./ sum(CN);
iem_noise = sum(CN) * tbin ./ max(CN);
iem99_noise = sum(CN) * tbin ./ prctile(CN, percentile);
plot([tm_noise', iem_noise', iem99_noise'])
legend({'TM', 'IEM', 'IEM99'})
xlabel('Sample');
ylabel('Lifetimes (ps)')

% SNR
S_tm = mean(tm_noise(l+1:end)) - mean(tm_noise(1:l));
S_iem = mean(iem_noise(l+1:end)) - mean(iem_noise(1:l));
S_iem99 = mean(iem99_noise(l+1:end)) - mean(iem99_noise(1:l));
N_tm = mean([std(tm_noise(l+1:end)), std(tm_noise(1:l))]);
N_iem = mean([std(iem_noise(l+1:end)), std(iem_noise(1:l))]);
N_iem99 = mean([std(iem99_noise(l+1:end)), std(iem99_noise(1:l))]);
title(sprintf('Noise = %0.1f. SNR: Tm-%0.2f, IEM-%0.2f, IEM%i-%0.2f\n', sigma,...
    S_tm/N_tm, S_iem/N_iem, percentile, S_iem99/N_iem99));

%% Test different noise levels
% Conclusion: IEM99 has a bit higher SNR than IEM at all noise levels. Tm
% is the the king of SNR at all levels. But, at very high sigma value, the
% difference is diminishing between the 3 metrics.

% Experimental parameters
dff = 2;
l = 100;

% IEM99
percentile = 99;

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;
C = a2 * i1 + (b2 * i2) * dff;

% Generate noise
sigmas = 0.05 : 0.05 : 0.7;
SNRs = zeros(length(sigmas), 3);

for i = 1 : length(sigmas)
    % Noise
    sigma = sigmas(i);
    N = randn(size(C)) * sigma + 1;
    N(N < 0) = 0;
    CN = C .* N;
    
    % New lifetimes
    tm_noise = sum(CN .* (t * ones(1, l*2)), 1) ./ sum(CN);
    iem_noise = sum(CN) * tbin ./ max(CN);
    iem99_noise = sum(CN) * tbin ./ prctile(CN, percentile);
    
    
    % SNR
    S_tm = mean(tm_noise(l+1:end)) - mean(tm_noise(1:l));
    S_iem = mean(iem_noise(l+1:end)) - mean(iem_noise(1:l));
    S_iem99 = mean(iem99_noise(l+1:end)) - mean(iem99_noise(1:l));
    N_tm = mean([std(tm_noise(l+1:end)), std(tm_noise(1:l))]);
    N_iem = mean([std(iem_noise(l+1:end)), std(iem_noise(1:l))]);
    N_iem99 = mean([std(iem99_noise(l+1:end)), std(iem99_noise(1:l))]);
    SNRs(i,:) = [S_tm/N_tm, S_iem/N_iem, S_iem99/N_iem99];
end

% Plot
figure
plot(sigmas', SNRs);
xlabel('Noise sigma');
ylabel('SNR');
legend({'TM', 'IEM', 'IEM99'})
title('Noise sigmas')

%% Test different percentiles
% Conclusion: IEM99 should probably be the 98-99th percentile. Anything
% less is just adding noise

% Experimental parameters
dff = 2;
l = 100;

% IEM99
percentiles = 94 : 0.5 : 100;

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;
C = a2 * i1 + (b2 * i2) * dff;

% Generate noise
sigma = 0.3;
SNRs = zeros(length(percentiles), 3);

for i = 1 : length(percentiles)
    % Percentile
    percentile = percentiles(i);
    
    % Noise
    N = randn(size(C)) * sigma + 1;
    N(N < 0) = 0;
    CN = C .* N;
    
    % New lifetimes
    tm_noise = sum(CN .* (t * ones(1, l*2)), 1) ./ sum(CN);
    iem_noise = sum(CN) * tbin ./ max(CN);
    iem99_noise = sum(CN) * tbin ./ prctile(CN, percentile);
    
    
    % SNR
    S_tm = mean(tm_noise(l+1:end)) - mean(tm_noise(1:l));
    S_iem = mean(iem_noise(l+1:end)) - mean(iem_noise(1:l));
    S_iem99 = mean(iem99_noise(l+1:end)) - mean(iem99_noise(1:l));
    N_tm = mean([std(tm_noise(l+1:end)), std(tm_noise(1:l))]);
    N_iem = mean([std(iem_noise(l+1:end)), std(iem_noise(1:l))]);
    N_iem99 = mean([std(iem99_noise(l+1:end)), std(iem99_noise(1:l))]);
    SNRs(i,:) = [S_tm/N_tm, S_iem/N_iem, S_iem99/N_iem99];
end

% Plot
figure
plot(percentiles', SNRs);
xlabel('Percentiles');
ylabel('SNR');
legend({'TM', 'IEM', 'IEM99'})
title(sprintf('Noise sigma = %0.1f. Percentiles for IEM99', sigma))

%% Test different bin sizes
% Conclusion: When bin sizes are small. IEM99 and Tm become really good at
% SNR. For IEM, it's probably because the 99th percentile estimation
% becomes very accurate. Where as the max value for IEM is noisy no matter
% what. Run the first block again after this since the bin sizes have been
% changed.

% Bins
tbins = cat(2, 1:2:8, 10 : 10 : 100); % In ps

% Experimental parameters
dff = 2;
l = 100;

% Percentile
percentile = 99;

% Noise
sigma = 0.3;
SNRs = zeros(length(tbins), 3);

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;

for i = 1 : length(tbins)
    % Get bin
    tbin = tbins(i);

    t = (1 : tbin: 10000)'; % in units of ps
    atau = 500; % asum is also tau_a, in ps
    btau = 2000; % sum is also tau_b, in ps
    a = exp(-t/atau); % Lifetime distribution of a
    b = exp(-t/btau); % Lifetime distribution of b
    fprintf('Tbin = %i. Sum a-%0.1f b-%0.1f\n', tbin, sum(a), sum(b));
    
    % Normalize it out
    asum = sum(a);
    bsum = sum(b);
    a2 = a/asum*100;
    b2 = b/bsum*100;
    
    % Generate trace
    C = a2 * i1 + (b2 * i2) * dff;
    
    % Noise
    N = randn(size(C)) * sigma + 1;
    N(N < 0) = 0;
    CN = C .* N;

    % New lifetimes
    tm_noise = sum(CN .* (t * ones(1, l*2)), 1) ./ sum(CN);
    iem_noise = sum(CN) * tbin ./ max(CN);
    iem99_noise = sum(CN) * tbin ./ prctile(CN, percentile);
    
    % SNR
    S_tm = mean(tm_noise(l+1:end)) - mean(tm_noise(1:l));
    S_iem = mean(iem_noise(l+1:end)) - mean(iem_noise(1:l));
    S_iem99 = mean(iem99_noise(l+1:end)) - mean(iem99_noise(1:l));
    N_tm = mean([std(tm_noise(l+1:end)), std(tm_noise(1:l))]);
    N_iem = mean([std(iem_noise(l+1:end)), std(iem_noise(1:l))]);
    N_iem99 = mean([std(iem99_noise(l+1:end)), std(iem99_noise(1:l))]);
    SNRs(i,:) = [S_tm/N_tm, S_iem/N_iem, S_iem99/N_iem99];
end

% Plot
figure
plot(tbins', SNRs);
xlabel('Bin sizes (ps)');
ylabel('SNR');
legend({'TM', 'IEM', 'IEM99'})
title(sprintf('Noise sigma = %0.1f. Testing bin sizes.', sigma))

%% Test photon number contamination (e.g., freaky laser, not due to dff)
% Conclusion: no obvious impact.

% Freaky laser power fluctuations
lasers = 1 : 0.5 : 10;

% Experimental parameters
dff = 2;
l = 100;

% fractions
i1 = cat(2, ones(1,l)*10, ones(1,l));
i2 = 11 - i1;
C = a2 * i1 + (b2 * i2) * dff;

% Percentile
percentile = 99;

% Noise
sigma = 0.3;
SNRs = zeros(length(lasers), 3);

for i = 1 : length(lasers)
    % Laser
    laser = lasers(i);

    % Noise
    N = randn(size(C)) * sigma + 1;
    N(N < 0) = 0;
    CN = C * laser .* N;

    % New lifetimes
    tm_noise = sum(CN .* (t * ones(1, l*2)), 1) ./ sum(CN);
    iem_noise = sum(CN) * tbin ./ max(CN);
    iem99_noise = sum(CN) * tbin ./ prctile(CN, percentile);
    
    % SNR
    S_tm = mean(tm_noise(l+1:end)) - mean(tm_noise(1:l));
    S_iem = mean(iem_noise(l+1:end)) - mean(iem_noise(1:l));
    S_iem99 = mean(iem99_noise(l+1:end)) - mean(iem99_noise(1:l));
    N_tm = mean([std(tm_noise(l+1:end)), std(tm_noise(1:l))]);
    N_iem = mean([std(iem_noise(l+1:end)), std(iem_noise(1:l))]);
    N_iem99 = mean([std(iem99_noise(l+1:end)), std(iem99_noise(1:l))]);
    SNRs(i,:) = [S_tm/N_tm, S_iem/N_iem, S_iem99/N_iem99];
end

% Plot
figure
plot(lasers', SNRs);
xlabel('Laser power multiplier');
ylabel('SNR');
legend({'TM', 'IEM', 'IEM99'})
title(sprintf('Noise isgma = %0.1f. Testing laser powers.', sigma))
