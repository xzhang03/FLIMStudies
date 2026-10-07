%% FLIM study 3 new
% Stephen X. Zhang (2023/06/14)
% The goal is to understand how fluorescence intensity changes may affect
% lifetime estimates. This is because adjusting amplitudes will have
% non-linear scaling effects on intensities, depending on tau (i.e., *the
% amplitude-intensity scaling problem*). This simulation treats intensity
% and lifetime as two independent properties.

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

% Updated 2024/9/19

clear
t = (1 : 10000)'; % in units of ps
asum = 500; % asum is also tau_a, in ps
bsum = 2000; % sum is also tau_b, in ps
a = exp(-t/asum); % Lifetime distribution of a
b = exp(-t/bsum); % Lifetime distribution of b

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

%% Add IRF (optional)
% Not necessary for this simulation but you can try it.
% IRF is assumed to be a simple gaussian
%{
sigma = 30;
irfx = (-3 * sigma : 3 * sigma)';
irf = exp(-(irfx / sigma) .^ 2 / 2);

a = conv(a,irf) / sum(irf);
b = conv(b,irf) / sum(irf);
t = (1 : 10000 + length(irf) - 1)';
%}

%% Situation 1: A-B transtiion, both 100% emission (FLIM akar)
% This scenario tests when the two states are qually bright 
% In equal-molecule condition: intensity a = intensty b
% The main conclusions: 
% 1. Tm scales linearly with the fraction of b.
% 2. IEM scales supralinearly with the fraction of b. 
% 3. If we don't have knowlwdge of fraction b, and just plot Tm against IEM,
% Tm will plateau at high IEM (high b)

% Testing a mixture of a and b with varying ratios
% i1 = amplitude of a
% i2 = amplitude of b
% i1 + i2 not = 1 due to *the amplitude-intensity scaling problem*

i1 = 1 : -0.1 : 0.1; % Varying amplitude a
i2 = (0.1 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
C = a * i1 + b * i2; % The combine distritbuion of a + b (detected at PMT)
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% IEM and Tm
tm = sum(C .* (t * ones(1,10))) ./ sum(C);
iem = sum(C) ./ max(C, [], 1);

figure('Name', 'Situation 1: a and b equal brightness')
subplot(1,3,1)
plot(fb, tm, '-o')
xlabel('Fraction of B')
title('Tm')
ylims = get(gca, 'YLim');

subplot(1,3,2)
plot(fb, iem, '-o')
xlabel('Fraction of B')
title('IEM')
ylim(ylims);

subplot(1,3,3)
plot(iem, tm, '-o')
xlabel('IEM')
ylabel('Tm')
title('IEM vs Tm')
xlim(ylims);
ylim(ylims);

%% Situation 2: A-B transition, short tau species 25% emission (long-bright: dlight/rcamp)
% This scenario tests when b state (long tau) is 4x brighter than a state
% (short tau). dff = 4
% In equal-molecule condition: intensity a = intensty b / dff
% The main conclusions: 
% 1. Tm scales sublinearly with the fraction of b.
% 2. IEM scales linearly with the fraction of b. 
% 3. If we don't have knowlwdge of fraction b, and just plot Tm against IEM,
% Tm will plateau at high IEM (high b)

% Testing a mixture of a and b with varying ratios
% i1 = amplitude of a
% i2 = amplitude of b
% i1 + i2 not = 1 due to *the amplitude-intensity scaling problem*
% To simulate dff, we will assume every b molecule contributes photon,
% while 1 out of 4 a molecules contributes photon (the protein is still there just
% dark)

dff = 4; % Defined as Fb/Fa (ideal measurements). Calling it dff is a bit wrong here since it doesn't subtract 1.

i1 = 1 : -0.1 : 0.1; % Varying amplitude a
i2 = (0.1 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
C = (a * i1) + (b * i2) * dff; % The combine distritbuion of a + b (detected at PMT), taking dff into account
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% IEM and Tm
tm = sum(C .* (t * ones(1,10))) ./ sum(C);
iem = sum(C) ./ max(C, [], 1);

figure('Name', 'Situation 2: b (long tau form) 4x brighter')

subplot(1,3,1)
plot(fb, tm, '-o')
xlabel('Fraction of B')
title('Tm')
ylabel('Tm (ps)')

subplot(1,3,2)
plot(fb, iem, '-o')
xlabel('Fraction of B')
title('IEM')
ylabel('IEM (ps)')

ylims = get(gca, 'YLim');
subplot(1,3,1)
ylim(ylims);

subplot(1,3,3)
plot(iem, tm, '-o')
xlabel('IEM (ps)')
ylabel('Tm (ps)')
title('IEM vs Tm')
xlim(ylims);
ylim(ylims);

%% Situation 3: a-b transition, long tau species 25% emission (short-bright).
% This scenario tests when b state (long tau) is 4x brighter than a state
% (short tau). dff = 0.25. I can't think of sensors with this property but
% testing it just in case.
% In equal-molecule condition: intensity a = intensty b / dff
% The main conclusions: 
% 1. Tm scales supralinearly with the fraction of b.
% 2. IEM scales even more supralinearly with the fraction of b and has
% small dynamic range
% 3. If we don't have knowlwdge of fraction b, and just plot Tm against IEM,
% Tm has better dynamic range than IEM

% Testing a mixture of a and b with varying ratios
% i1 = amplitude of a
% i2 = amplitude of b
% i1 + i2 not = 1 due to *the amplitude-intensity scaling problem*
% To simulate dff, we will assume every a molecule contributes photon,
% while 1 out of 4 b molecules contributes photon (the protein is still there just
% dark)

dff = 0.25; % Defined as Fb/Fa (ideal measurements). Calling it dff is a bit wrong here since it doesn't subtract 1.

i1 = 1 : -0.1 : 0.1; % Varying amplitude a
i2 = (0.1 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
C = (a * i1) + (b * i2) * dff; % The combine distritbuion of a + b (detected at PMT), taking dff into account
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% IEM and Tm
tm = sum(C .* (t * ones(1,10))) ./ sum(C);
iem = sum(C) ./ max(C, [], 1);

figure('Name', 'Situation 3: a (short tau form) 4x brighter')

subplot(1,3,1)
plot(fb, tm, '-o')
xlabel('Fraction of B')
title('Tm')
ylims = get(gca, 'YLim');

subplot(1,3,2)
plot(fb, iem, '-o')
xlabel('Fraction of B')
title('IEM')
ylim(ylims);

subplot(1,3,3)
plot(iem, tm, '-o')
xlabel('IEM')
ylabel('Tm')
title('IEM vs Tm')
xlim(ylims);
ylim(ylims);

%% Demonstration 1. Tm dose-response curves with different intensity changes (B>A)
% Figure
figure('Name', 'Demonstration 1:Tm dose-response curve with different intensity changes (B>A)')

% New plot
new = true;

% Series
i1 = 1 : -0.1 : 0; % Varying amplitude a
i2 = (0 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% dffs
dffs = [1 2 4 8 16];

% Defined as Fb/Fa (ideal measurements). Calling it dff is a bit wrong here since it doesn't subtract 1.
for dff = dffs
    % Plot
    C = (a * i1) + (b * i2) * dff; % The combine distritbuion of a + b (detected at PMT), taking dff into account
    
    % IEM and Tm
    tm = sum(C .* (t * ones(1,11))) ./ sum(C);
    
    hold on
    plot(fb, tm, '-o')
    hold  off
    
    if new 
        new = false;
        xlabel('Fraction of B (X)')
        title('Intensity-weighted averaging (Tm), B brighter than A')
        ylabel('Tm (ps)')
    end
end

legends = cell(1,length(dffs));
for i = 1 : length(dffs)
    legends{i} = sprintf('dF/F0 = %i', dffs(i)- 1);
end

legend(legends)

%% Demonstration 2. Tm dose-response curves with different intensity changes (B<A, rare)
% Figure
figure('Name', 'Demonstration 2: Tm dose-response curves with different intensity changes (B<A, rare)')

% New plot
new = true;

% Series
i1 = 1 : -0.1 : 0; % Varying amplitude a
i2 = (0 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% dffs
dffs = [1 0.5 0.25 0.125 0.0625];

% Defined as Fb/Fa (ideal measurements). Calling it dff is a bit wrong here since it doesn't subtract 1.
for dff = dffs
    % Plot
    C = (a * i1) + (b * i2) * dff; % The combine distritbuion of a + b (detected at PMT), taking dff into account
    
    % IEM and Tm
    tm = sum(C .* (t * ones(1,11))) ./ sum(C);
    
    hold on
    plot(fb, tm, '-o')
    hold  off
    
    if new 
        new = false;
        xlabel('Fraction of B (X)')
        title('Intensity-weighted averaging (Tm), B dimmer than A (very rare)')
        ylabel('Tm (ps)')
    end
end

legends = cell(1,length(dffs));
for i = 1 : length(dffs)
    legends{i} = sprintf('dF/F0 = -%0.2f', 1-dffs(i));
end

legend(legends)

%% Demonstration 3. IEM targeted simulation
% Figure
figure('Name', 'Demonstration 3: IEM targeted simulation')

% New plot
new = true;

% Series
i1 = 1 : -0.1 : 0; % Varying amplitude a
i2 = (0 : 0.1 : 1) / bsum * asum; % Calculating the amplitude b such that molecule-counts don't change
fb = (i2 * bsum) ./ (i1 * asum + i2 * bsum); % Calculating the fraction of B for plotting purposes

% dffs
dffs = [1 2 4 8 16];

% Defined as Fb/Fa (ideal measurements). Calling it dff is a bit wrong here since it doesn't subtract 1.
for dff = dffs
    % Plot
    C = (a * i1) + (b * i2) * dff; % The combine distritbuion of a + b (detected at PMT), taking dff into account
    
    % IEM and Tm
    iem = sum(C) ./ max(C, [], 1);
    
    hold on
    plot(fb, iem, '-o')
    hold  off
    
    if new 
        new = false;
        xlabel('Fraction of B (X)')
        title('Amplitude-weighted averaging (IEM), TauB/TauA = 4')
        ylabel('IEM (ps)')
    end
end

legends = cell(1,length(dffs));
for i = 1 : length(dffs)
    legends{i} = sprintf('Ir(B/A) = %i', dffs(i));
end

legend(legends)

