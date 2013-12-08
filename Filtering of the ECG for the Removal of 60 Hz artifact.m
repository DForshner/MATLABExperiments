% Filtering of the ECG for the Removal of 60 Hz artifact
%
% Develop a MATLAB program to implement the notch filter designed as above. Apply
% the filter to the signal in the data file ecg 60hz 200.dat.
%
% Apply the notch filter to a segment of one of your noisy ECG signals from Lab 1 with
% 3−5 cardiac cycles. Study the nature of the artifacts in the noisy signal and in the output
% of the filter.
%
% Analyze the characteristics and the effects of the filter in the frequency domain by
% obtaining the pole-zero diagram and the transfer function of the filter, as well as the Fourier
% spectra of the input and output signals.

close all; clear all; clc; % Clear everything

% -------------------------------------------------------------- Part 1
% Read and plot the ECG

ECG_Voltage = load('ecg_60hz_200.dat'); % Load ECG data

fs = 200; % Sampling frequency is fs = 200 Hz
N = length(ECG_Voltage); % number of samples
T = 1/fs; % period
ECG_Time = [0 : (N - 1)] * T; % construct time axis

%Plot the ECG Signal vs. Time
ECG1Plot = figure('Name','ecg60hz200.dat - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECG_Time, ECG_Voltage);
title('ecg60hz200.dat - ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (Unknown)');
axis tight;

% Determine the where the zeros go on the notch filter
OmegaZero = 2 * pi * (60 / 200) % OmegaZero = 2*pi (notch freq / sampling freq)
z1 = cos(OmegaZero) + 1i * sin(OmegaZero)
z2 = cos(-OmegaZero) + 1i * sin(-OmegaZero)

% Convert the transfer function H(z) to a difference equation y(n)
% H(z) = (z-z1)(z-z2) / z^2 
% H(z) = (Z^2 - (z1+z2)z + z1z2) / z^2
% H(z) = 1 - (z1+z2)z^-1 + (z1z2)z^-2
% H(z) = Y(z) / X(z)
% Y(z) = X(z) * (1 - (z1+z2)z^-1 + (z1z2)z^-2)
% Y(z) = X(z) - X(z)(z1+z2)z^-1 + X(z)(z1z2)z^-2
% y(n) = x(n) - (z1+z2)*x(n-1) + (z1z2)*x(n-2)
% y(n) = b(1)x(n) + b(2)x(n-1) + b(3)x(n-2)

a = []; b = [];

a(1) = 1
b(1) = 1;
b(2) = -z1 + -z2;
b(3) = z1 * z2

% For low-pass and notch filters we set normzlize the DC gain to 1
% H(z) = Gain * (b(1) + b(2)z^-1 + b(3)z^-2)
% H(z=1) = Gain * (b(1) + b(2)1 + b(3)1) = 1
HZEqualOne = b(1) + b(2) + b(3);
Gain = 1 / HZEqualOne
b = b * Gain

% Filter the ECG signal
ECGFiltered_Voltage = filter(b, a, ECG_Voltage);

%Plot the Filtered ECG Signal vs. Time
subplot(2,1,2); plot(ECG_Time, ECGFiltered_Voltage);
title('ecg60hz200.dat - Filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (Unknown)');
axis tight;

% Plot the PSD of the unfiltered signal
ECG1PSGPlot = figure('Name','ecg60hz200.dat - PSDs of ECG'); % Create a new figure
subplot(2,1,1);
ECG = FindPSD2(ECG_Voltage, (length(ECG_Voltage)), fs);
title('ecg60hz200.dat - PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (fs/2) -40 5]); % If leave to auto it will be hard to compare

% Plot the PSD of the filtered signal
subplot(2,1,2);
ECGFiltered = FindPSD2(ECGFiltered_Voltage, (length(ECGFiltered_Voltage)), fs);
title('ecg60hz200.dat - PSD of filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (fs/2) -40 5]); % If leave to auto it will be hard to compare

% -------------------------------------------------------------- Part 2
% Read and plot the ECG
load SampleECG.txt;  % Load file into memory
ECG_Time = SampleECG(1:520,1);   %  copy column 1 of the matrix to Time
ECG_Voltage = SampleECG(1:520,2);   %  copy column 2 of the matrix to Voltage

%Plot the ECG Signal vs. Time
ECG2Plot = figure('Name','SampleECG.txt - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECG_Time, ECG_Voltage);
title('SampleECG.txt - ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (Unknown)');
axis tight;

% Filter the ECG signal
ECG2Filtered_Voltage = filter(b, a, ECG_Voltage);

%Plot the Filtered ECG Signal vs. Time
subplot(2,1,2); plot(ECG_Time, ECG2Filtered_Voltage);
title('SampleECG.txt - Filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (Unknown)');
axis tight;

% Plot the PSD of the unfiltered signal
ECG1PSGPlot = figure('Name','SampleECG.txt - PSDs of ECG'); % Create a new figure
subplot(2,1,1);
ECG2 = FindPSD2(ECG_Voltage, (length(ECG_Voltage)), fs);
title('SampleECG.txt - PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (fs/2) -40 5]); % If leave to auto it will be hard to compare

% Plot the PSD of the filtered signal
subplot(2,1,2);
ECG2Filtered = FindPSD2(ECG2Filtered_Voltage, (length(ECG2Filtered_Voltage)), fs);
title('SampleECG.txt - PSD of Filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (fs/2) -40 5]); % If leave to auto it will be hard to compare

% -------------------------------------------------------------- Part 3
% Plot the filter responce
FilterFreqResponse = figure('Name','Filter Frequency & Phase Response'); % Create a new figure
freqz(b,a, 512, 200);

% Plot the filter pole - zero diagram
ZPlaneFilter = figure('Name','Pole-Zero Diagram of the transfer function of the filter'); % Create a new figure
freqz(b,a, 512, 200);
zplane(b,a)