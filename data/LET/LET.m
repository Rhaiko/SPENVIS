% Clear workspace.
clc;
close all
clear;

% Set format as long so no data is lost when reading the data.
format long

% Read the differential and integrated flux table from SPENVIS.
L = importdata('LET.txt');

% Make data readable.
x = str2double(L.rowheaders);
y1 = L.data(:,1);
y2 = L.data(:,2);

% Integrate the flux by hand.
y3 = zeros(length(y1),1);
for i = 1:length(y2) -1
    y3(i + 1) = sum(y2(length(y2) - i:end));
end

% Manipulate (flip) data.
y4 = flipud(y3);
y5 = (4 * pi) .* y4;

% Plot data.
figure;
loglog(x(2:end), y2(2:end), x(2:end), y1(2:end), x(2:end), y5(2:end))
legend({'Differential flux SPENVIS [(m^{-2} sr^{-1} s^{-1})(MeV cm^2 g^{-1})^{-1}]', 'Integral flux SPENVIS [m^{-2} sr^{-1} s^{-1}]', 'Integral flux [m^{-2} sr^{-1} s^{-1}]'}, 'Location', 'southwest')
xlabel('LET(Si) [MeV cm^2 g^{-1}]')
