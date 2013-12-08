% Synchronized Averaging for Noise Reduction
% 
% We have been provided with M sets of ERP signals each with N data points
% sampled at a rate of 1000 Hz.  This code Plots the ERP signals for the 
% first two cases(E11 to E44 and E11 to E88) and the synchronized average 
% signal for all cases.  Also calculates the
% SNR & Euclidean Distance for each case. (Portions of the code are
% directly copied from esopX.m)

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

% Read signal data in the order of stimulation into 2 dim vectors
load E11; load E22; load E33; load E44; load E55; load E66
load E77; load E88; load E99; load E1010; load E1111
load E1212; load E1313; load E1414; load E1515; load E1616
load E1717; load E1818; load E1919; load E2020; load E2121
load E2222; load E2323; load E2424

% Copy the ERP data vectors into 3 dim data vector d
d(:,1) = E11; d(:,2) = E22; d(:,3) = E33; d(:,4) = E44; d(:,5) = E55;
d(:,6) = E66; d(:,7) = E77; d(:,8) = E88; d(:,9) = E99; d(:,10) = E1010;
d(:,11) = E1111; d(:,12) = E1212; d(:,13) = E1313; d(:,14) = E1414;
d(:,15) = E1515; d(:,16) = E1616; d(:,17) = E1717; d(:,18) = E1818;
d(:,19) = E1919; d(:,20) = E2020; d(:,21) = E2121; d(:,22) = E2222;
d(:,23) = E2323; d(:,24) = E2424;

% Setup the parameters of the data acquisition
fs = 1000;     % sampling rate
gain = 10000;  % amplifier gain
T = 1.0/fs;    % sampling interval in seconds
NMax = size(d,1); % number of samples in each signal
t = [1:NMax] * T; % time axis in seconds: common to all signals

% ------------------------------------------------------- E11 to E44
%
% Create a new figure and plot the signals and their sync average
figure_E11toE44 = figure('Name','Signals E11toE44 and their synchronized Average');
disp('E11toE44')

% Setup the samples to be averaged
k_min = 1;              % first signal to be plotted
k_max = 4;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

plot_number = 1;
for k = k_min:k_max

    % Plot the current signal in the proper subplot
    subplot((MMax + 1), 1, plot_number);
    plot(t, d(:,k), 'k-');
    Title(['Signal m= ', num2str(k)]);
    xlabel('Time (sec)');
    ylabel('ERP*10^4V');
    axis tight;  
    plot_number = plot_number + 1;
end

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot((MMax + 1), 1, (MMax + 1));
plot(t, yofn,'k-');
TITLE('Sync Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt(sum((d(:,M) - yofn).^2));   
end
EuclideanDistance = runningM_sum / MMax

% ------------------------------------------------------- E11 to E88
%
% Create a new figure and plot the signals and their synchronized average
disp('E11toE88')
figure_E11toE88 = figure('Name','Signals E11toE88 and their synchronized Average');

% Setup the samples to be averaged
k_min = 1;              % first signal to be plotted
k_max = 8;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

plot_number = 1;
for k = k_min:k_max
    % Plot the current signal in the proper subplot
    subplot( (MMax + 1), 1, plot_number);
    plot(t, d(:,k), 'k-');
    Title(['Signal m= ', num2str(k)]);
    xlabel('Time (sec)');
    ylabel('ERP*10^4V');
    axis tight;  
    plot_number = plot_number + 1;
end

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot((MMax + 1), 1, (MMax + 1));
plot(t, yofn,'k-');
TITLE('Sync Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max  
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt( sum( (d(:,M) - yofn).^2 ) );   
end
EuclideanDistance = runningM_sum / MMax

% Prepare a new figure that will plot the 4 remaining Sychronized Averages
figure_Add_Sync_Ave = figure('Name','Additional Synchronized Averages');

% ------------------------------------------------------- E11 to E1212
%
disp('E11toE1212')

% Setup the samples to be averaged
k_min = 1;              % first signal to be plotted
k_max = 12;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot(4, 1, 1);
plot(t, yofn,'k-');
TITLE('E11toE1212 - Synchronized Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt( sum( (d(:,M) - yofn).^2 ) );   
end
EuclideanDistance = runningM_sum / MMax

% ------------------------------------------------------- E11 to E2424
%
disp('E11toE2424')

% Setup the samples to be averaged
k_min = 1;              % first signal to be plotted
k_max = 24;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot(4, 1, 2);
plot(t, yofn,'k-');
TITLE('E11toE2424 - Synchronized Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt( sum( (d(:,M) - yofn).^2 ) );   
end
EuclideanDistance = runningM_sum / MMax

% ------------------------------------------------------- E1717 to E2424
%
disp('E1717 to E2424')

% Setup the samples to be averaged
k_min = 17;              % first signal to be plotted
k_max = 24;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot(4, 1, 3);
plot(t, yofn,'k-');
TITLE('E1717 to E2424 - Synchronized Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max   
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt( sum( (d(:,M) - yofn).^2 ) );   
end
EuclideanDistance = runningM_sum / MMax

% ------------------------------------------------------- E1313 to E2424
%
disp('E1313 to E2424')

% Setup the samples to be averaged
k_min = 13;              % first signal to be plotted
k_max = 24;              % last signal to be plotted
MMax = k_max - k_min + 1;   % number of signals to be plotted
disp(['M = ' , num2str(MMax)]); % Display the number of samples in the console

% Calculate and plot the sychronized average signal
yofn=sum(d(:,k_min:k_max),2)/MMax;
subplot(4, 1, 4);
plot(t, yofn,'k-');
TITLE('E1313 to E2424 - Synchronized Average');
xlabel('Time(sec)');
ylabel('ERP*10^4V');
axis tight;

% Calculate Noise Power
runningM_sum = 0;
for M = k_min:k_max   
    runningM_sum = runningM_sum + sum((d(:,M) - yofn).^2);
end
NoisePowerSquared = runningM_sum / (NMax * 0.001 *(MMax - 1));

% Calculate Signal Power
SignalPowerSquared = sum(yofn.^2)/(NMax * 0.001) - (NoisePowerSquared / MMax);

% Calculate Signal to Noise Ratio
SNR = SignalPowerSquared / NoisePowerSquared

% Calculate Euclidean Distance
runningM_sum = 0;
for M = k_min:k_max
    runningM_sum = runningM_sum + sqrt( sum( (d(:,M) - yofn).^2 ) );   
end
EuclideanDistance = runningM_sum / MMax