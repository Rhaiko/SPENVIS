% Clear workspace.
clear
close all
clc

% The mono-beam experimental results of SMJ329C50GFAM66:
% E - Energy of Beam    - MeV
% F - Flux of particles - cm^-2 s^-1
% t - Exposure time     - s
% n - Number of flips   - Arbitrary
% Ions: C12    C12    C12    O16   Ar40   Fe56   Kr84  Xe131
E   = [ 0.6,  0.72,   9.6,   4.8,    20,    56,    84,   786];
F   = [  25,    25,    25,    17,    23,    10,     7,     4].*(10^6);
t   = [   5,     5,     5,     5,     5,    10,    10,    15].*60;
n   = [   8,  7497, 22514, 23986, 33810, 29991, 21022, 18043];

% Extracted from CREME stopping powers in Silicon [MeV cm^2 mg^-1].
LET = [ 2.5,   3,   4.5,     7,      18,    29,    40,   60];

% Cross section.
sig = n./(F.*t); % cm^2

% Plot sigma to LET.
figure;
plot(LET, sig)
xlabel('LET [MeV cm^2 mg^-1]')
ylabel('sigma [m^2]')
% L0 approx 4.5 MeV cm^2 mg^-1 (threshold in LET).
% Cs approx 5 * 10^-6 m^2 (saturated cross section).

% Set format as long so no data is lost when reading the data.
format long

% Read the differential and integrated flux table from SPENVIS.
L = importdata('LET.txt');
flux(:,1) = str2double(L.rowheaders);
flux(:,2) = L.data(:,1);
flux(:,3) = L.data(:,1);

% Cruve fitting.
weibull = @(v,x) v(1).*(1-exp(-((x-v(2))./v(3)).^v(4)));
options = optimset('TolFun', 1e-60);
v       = abs(lsqcurvefit(weibull, [max(sig), 4.5, 5.7, 1], LET, sig));

sigma = zeros(1000,4);
sigma(:,4) = (flux(:,1)>v(2)).*v(1).*(1-exp(-((flux(:,1)-v(2))/v(3)).^v(4)));

% chip:  nmos       cmos  bipolar
param = [1.71e-5, 1.2e-5,  2.6e-5;      % L0 - MeV cm^2 mg^-1
           0.487,  136.8,     0.6;      % Cs - cm^2
            4.95,    350,     4.4;      % W  - MeV cm^2 mg^-1
           1.422,      3,     0.7];     % s  - arbirtary

% Calculate sigma at each moment.
for i = 1:3
    sigma(:,i) = (flux(:,1)>param(2,i)).*param(1,i).*(1-exp(-((flux(:,1)-param(2,i))/(param(3,i))).^param(4,i)));
end

% Plot sensitivity.
figure;
for i = 1:4
    loglog(flux(:,1),sigma(:,i))
    hold on
end
title('Sensitivity')
ylabel('sigma(LET) [cm^2]')
xlabel('LET [MeV cm^2 mg^{-1}]')
xlim([0 1.5*10^3])
ylim([10^-7 10^-4])
legend('NMOS', 'CMOS', 'Bipolar', 'SMJ', 'Location', 'southeast')
hold off

% Single Event Upset Calculation per chip.
% NMOS; CMOS; Bipolar; SMJ.
SEU = zeros(4, 1);
for i = 1:4
    SEU(i) = 4*pi*trapz(flux(2:end,1), sigma(2:end,i).*flux(2:end,3))/999;
end








