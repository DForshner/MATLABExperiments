% Frequency-domain Analysis of Biomedical Signals
% Experiment 2
%
% Download one of the eeg1*.dat files; there are eight channels available. The sampling rate is
% fs = 100 Hz. Obtain and analyze the PSD of the entire EEG signal in the selected channel.
% Include in your report a plot of the EEG signal selected and the corresponding PSD.

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

% ---------------------------------------------------- Experiment 2
% download one of the eeg1*.dat files; there are eight channels available. 
% The sampling rate is fs = 100 Hz. Obtain and analyze the PSD of the entire 
% EEG signal in the selected channel. Include in your report a plot of the 
% EEG signal selected and the corresponding PSD.

EEGSignal = load('eeg1-o1.dat');  % Load file into memory
EEGTimeAxis = (1 : length(EEGSignal)) * (1 / 100); % Sampling frequency is 100 Hz

EEGFigure = figure('Name','EEG signal of "eeg1-o1.dat" '); % Create a new figure

subplot(2,1,1), plot(EEGTimeAxis, EEGSignal)
title('EEG signal of "eeg1-o1.dat" ');
xlabel('Time (sec)');
ylabel('Amplitude (Unknown)');

EEG = FindPSD(EEGSignal, length(EEGSignal), 100);
subplot(2,1,2), hold on;
plot(EEG.HalfFs, EEG.PSD);
plot(EEG.HalfFs, -3, 'r'); % plot 3db line
hold off;
title('PSD of EEG signal of "eeg1-o1.dat" (fs = 100 Hz) (N = # samples in signal)');
xlabel('Frequency Hz');
ylabel('Magnitude (dB)');
axis([0 EEG.fs -40 5]); % If leave to auto it will be hard to compare