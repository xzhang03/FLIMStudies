%% FLIM study 4. Simulate TCSPC time jitter when estimating small lifetime changes
% The goal of this study is to determine the impact of TCSPC time jitter on
% measurements of small lifetime changes

% Conclusion: You NEED TCSPC time jitter to measure a small lifetime change
% (30 ps) when using large bins (50 ps). Noiseless TCSPCs are not good for
% this purpose. Jitter size needs to be >15 ps in the example above.

% Stephen Zhang, 7/19/2023

% Definition:
% Jitter: TCSPC imprecision in time measurements.
% Effect size: True lifetime changes
% Bin size: bin size of the TCSPC histogram

clear

% Suppose 30 ps signal, 2.5 ps jitter, and 50 ps bin-size
tm1 = 0;
tm2 = 30;
jitter = 2.5;
binsz = 50;

% Sampling
n_samples = 20000;
s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);

% Readout
(mean(s2) - mean(s1)) * binsz

%% Case 1. Jitter = 2.5 ps. Effectsize = 30 ps. Vary binsize from 1 to 100 ps
% Conclusion: Results are useless when binsize > 50
tm1 = 0;
tm2 = 30;
jitter = 2.5;
binszs = 1 : 2 : 100;

readouts = zeros(size(binszs));
for i = 1 : length(binszs)
    binsz = binszs(i);
    s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
    s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);
    readouts(i) = (mean(s2) - mean(s1)) * binsz;
end

plot(binszs, readouts)
xlabel('Bin size (ps)')
ylabel('Measured lifetime change (ps)')

%% Case 2. Jitter = 2.5 ps. binsize = 50 ps. Vary effect size from 10 to 200;
% Measured lifetime jumps in the following manner
% Step up at 25 (to 50), 75 (100), 125 (150
tm1 = 0;
jitter = 2.5;
tm2s = 10 : 10 : 200;
binsz = 50;

readouts = zeros(size(tm2s));
for i = 1 : length(tm2s)
    tm2 = tm2s(i);
    s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
    s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);
    readouts(i) = (mean(s2) - mean(s1)) * binsz;
end

plot(tm2s, readouts)
xlabel('Effect size (ps)')
ylabel('Measured lifetime change (ps)')

%% Case 3. Binsize = 50 ps. Effectsize = 30 ps. Vary jitter from 0 : 40
% Conclusion: To measure 30 ps difference using 50 ps bins, you need some 
% TCSPC jitter: about >15 ps. Jitter size is harmless in the 15-40 ps
% range.

tm1 = 0;
tm2 = 30;
jitters = 2 : 2 : 40;
binsz = 50;

readouts = zeros(size(jitters));
for i = 1 : length(jitters)
    jitter = jitters(i);
    s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
    s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);
    readouts(i) = (mean(s2) - mean(s1)) * binsz;
end

plot(jitters, readouts)
hold on
plot([jitters(1);jitters(end)], [tm2-tm1 tm2-tm1], '--k')
hold off
xlabel('TCSPC Jitter (ps)')
ylabel('Apparent lifetime (ps)')
title(sprintf('True lifetime = %i ps, Bin size = %i ps', tm2, binsz))

%% Case 4. 2D interaction between jitter size and bin size
% Conclusion: The only bad situation seems to be low jitter and large bin size
tm1 = 0;
tm2 = 30;
jitters = 4 : 4 :40;
binszs = 5 : 5 : 50;

readouts = zeros(length(jitters), length(binszs));
for i = 1 : length(jitters)
    jitter = jitters(i);
    for j = 1 : length(binszs)
        binsz = binszs(j);
        s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
        s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);
        readouts(i,j) = (mean(s2) - mean(s1)) * binsz - tm2;
    end
end

imagesc(readouts, [-20 20])

%% Case 5. Data following assay
% Mostly just a demonstration of the cyclical nature of errors in lifetime estimates.
tm1 = 0;
binsz = 50;
tm2s = 0 : 10 : 200;
jitters = 8 : 8 : 40;


readouts = zeros(length(jitters), length(tm2s));
for i = 1 : length(jitters)
    jitter = jitters(i);
    for j = 1 : length(tm2s)
        tm2 = tm2s(j);
        s1 = round((randn(n_samples,1) * jitter + tm1) / binsz);
        s2 = round((randn(n_samples,1) * jitter + tm2) / binsz);
        readouts(i,j) = (mean(s2) - mean(s1)) * binsz - tm2;
    end
end

plot(tm2s, readouts')
legend(cellfun(@num2str, num2cell(jitters), 'UniformOutput', false))
xlabel('True lifetime change (ps)')
ylabel('Error in estimating lifetime changes (ps)')