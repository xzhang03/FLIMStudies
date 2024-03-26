%% FLIM study 2. Simulate the photon-robustness and noise-robustness of FLIM analysis methods
% The goal of this simulation is to answer the question of how many photon
% is needed per histogram in noise-free and noise-added environment
% Stephen Zhang, 12/2/2023

% Simulate a sensor with one state: 2000 ps tau. The noise tau is 500 ps.

% Conclusions: 
% 1. TM needs >1.6k photons to be good (8x bin number) for absolute values. 
% 2. IEM needs >51k photons (256x bin number) for absolute values. This 
% makes it useful for bulk measurements such as 2p whole-FOV analysis and
% photometry. It's difficult to use for small ROIs

% Terminology
% Amplitude: the scaler factor before the expontial term
% Intensity: the area under curve 
% Tau: time constant of the exponential term
% [Intensity] = [Amplitude] * exp(t / tau);
% Tm: first moment ("mean arrival time")
% IEM: Ratio of Intensity/Amplitude

clear
tbin = 50;
t = (1 : tbin: 10000)'; % in units of ps
noisetau = 500; % asum is also tau_a, in ps
datatau = 2000; % sum is also tau_b, in ps
autofluo = round(exp(-t/noisetau)*10000); % Lifetime distribution of noise
data = round(exp(-t/datatau)*10000); % Lifetime distribution of b
fprintf('Sum noise-%0.1f data-%0.1f\n', sum(autofluo), sum(data));

% Find ideal values
datatm = sum(data .* t) / sum(data);
dataiem = sum(data) / max(data) * tbin;
fprintf('Data tm-%0.1f iem-%0.1f\n', datatm, dataiem);

figure
plot(t,autofluo,t,data);
xlabel('Time (ps)')
legend({'Autofluorescence (tau = 500 ps)'; 'Data (tau = 2000 ps)'})
title('Equal amplitude (not equal intensity)')

%% Case 1. Draw data of different sizes
% Conclusion: 
% 1. IEM needs at least 12k photons (>60x bin number) to converge mean,
% basically useful for absolute value estimates only in photometry/full FOV
% 2p experiments.
% 2. TM is reasonable above 1.6k photons (8x bin bin number) to converge mean.
% 3. RMS-wise, IEM needs 51k photons to be low 50 ps RMS.
% 4. TM is reaonsbale above 1.6 k photons in terms of RMS as well

% Drawing parameter
ns = 2.^(1:12)/4 * 200;
nrepeats = 20;

% Data matrix
tmmat = zeros(length(ns), nrepeats);
iemmat = zeros(length(ns), nrepeats);

for i = 1 : length(ns)
    % N
    n = ns(i);
    
    tic;
    for j = 1 : nrepeats
        % Use anonymouse function to draw data from distribution
        datadraw = drawhist(data, n);
        tm = mean(datadraw) * tbin;
        [~, pk] = mode(datadraw);
        iem = n / pk * tbin;

        % Save
        tmmat(i,j) = tm;
        iemmat(i,j) = iem;
    end
    
    % Time    
    ttoc = round(toc,1);
    if i == 1
        fprintf('Drawing samples: \n')
    end
    fprintf('%i: n = %i, T = %0.1f s\n', i, n, ttoc);
    if i == length(ns)
        disp('Done.')
    end
end

% Plot traces
figure
subplot(1,2,1)
plot(tmmat')
title('Tm traces')

subplot(1,2,2)
plot(iemmat')
title('IEM traces')

% Plot summary
figure

% Mean
subplot(1,2,1)
semilogx(ns, mean(tmmat,2))
hold on
semilogx(ns, mean(iemmat,2))
hold off
set(gca, 'XTick', ns)
title('Mean')
legend({'TM', 'IEM'})

% RMS
subplot(1,2,2)
semilogx(ns, rms(tmmat - mean(tmmat,2)*ones(1,nrepeats),2))
xlims = xlim();
hold on
semilogx(ns, rms(iemmat - mean(iemmat,2)*ones(1,nrepeats),2))
plot(xlims, [50 50], '--k')
hold off
set(gca, 'XTick', ns)
title('RMS')
legend({'TM', 'IEM'})

%% Case 2. Draw data of different sizes (added autofluorescence)
% Conclusion (same as above): 
% 1. IEM needs at least 12k photons (>60x bin number) to converge mean,
% basically useful for absolute value estimates only in photometry/full FOV
% 2p experiments.
% 2. TM is reasonable above 1.6k photons (8x bin bin number) to converge mean.
% 3. RMS-wise, IEM needs 51k photons to be low 50 ps RMS.
% 4. TM is reaonsbale above 1.6 k photons in terms of RMS as well

% Autofluorescence ratio
NR = 0.1;
data2 = round(data + autofluo * NR);

% Drawing parameter
ns = 2.^(1:12)/4 * 200;
nrepeats = 20;

% Data matrix
tmmat = zeros(length(ns), nrepeats);
iemmat = zeros(length(ns), nrepeats);

for i = 1 : length(ns)
    % N
    n = ns(i);
    
    tic;
    for j = 1 : nrepeats
        % Use anonymouse function to draw data from distribution
        datadraw = drawhist(data2, n);
        tm = mean(datadraw) * tbin;
        [~, pk] = mode(datadraw);
        iem = n / pk * tbin;

        % Save
        tmmat(i,j) = tm;
        iemmat(i,j) = iem;
    end
    
    % Time    
    ttoc = round(toc,1);
    if i == 1
        fprintf('Drawing samples: \n')
    end
    fprintf('%i: n = %i, T = %0.1f s\n', i, n, ttoc);
    if i == length(ns)
        disp('Done.')
    end
end

% Plot traces
figure
subplot(1,2,1)
plot(tmmat')
title('Tm traces')

subplot(1,2,2)
plot(iemmat')
title('IEM traces')

% Plot summary
figure

% Mean
subplot(1,2,1)
semilogx(ns, mean(tmmat,2))
hold on
semilogx(ns, mean(iemmat,2))
hold off
set(gca, 'XTick', ns)
title('Mean')
legend({'TM', 'IEM'})

% RMS
subplot(1,2,2)
semilogx(ns, rms(tmmat - mean(tmmat,2)*ones(1,nrepeats),2))
xlims = xlim();
hold on
semilogx(ns, rms(iemmat - mean(iemmat,2)*ones(1,nrepeats),2))
plot(xlims, [50 50], '--k')
hold off
set(gca, 'XTick', ns)
title('RMS')
legend({'TM', 'IEM'})


%% Anonymous functions
function hisout = drawhist(hisin, n)
    % Max index
    imax = sum(hisin);
    
    % cumulative his
    hiscdf = cumsum(hisin);

    % Draw
    rands = randi(imax, [n, 1]);

    % Generate output
    hisout = nan([n,1]);
    for i = 1 : n
        hisout(i) = find(rands(i) <= hiscdf, 1, 'first');
    end
end