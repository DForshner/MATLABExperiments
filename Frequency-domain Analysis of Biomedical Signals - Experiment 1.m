% Frequency-domain Analysis of Biomedical Signals
% Experiment 1
%
% (1) Select one cardiac cycle from your noise-free ECG recording from Lab 1. Obtain and analyze
% its PSD. Similarly, obtain the PSD of one cardiac cycle from your noisy
% ECG signal.
% (2) Repeat the experiment above with a larger number of samples in the FFT, such as 1024 or
% 2048, with zero padding of the signal;
% (3) Repeat the experiment with a segment containing 20 cycles of the ECG in your noise-free
% ECG recording;

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

% ---------------------------------------------------- Experiment 1 Part 1
% Obtain the PSDs of one cardiac cycle of clean and noisy ECG

OneCycleECGs=figure('Name','Clean and Noisy ECGs from Lab 1'); % Create a new figure

load CleanECG.txt;  % Load file into memory
CycleCleanECG_Time = CleanECG(310:455,1);   %  copy column 1 of the matrix to Time
CycleCleanECG_Voltage = CleanECG(310:455,2);   %  copy column 2 of the matrix to Voltage

subplot(2,1,1), plot(CycleCleanECG_Time,CycleCleanECG_Voltage)
title('ECG - 1 Complete Cycle of the clean signal');
xlabel('Time Frame (sec)'); ylabel('Voltage (mv)');

load MessyECG.txt;  % Load file into memory
CycleNoisyECG_Time = MessyECG(415:570,1);   %  copy column 1 of the matrix to Time
CycleNoisyECG_Voltage = MessyECG(415:570,2);   %  copy column 2 of the matrix to Voltage

subplot(2,1,2), plot(CycleNoisyECG_Time,CycleNoisyECG_Voltage)
title('ECG - 1 Complete Cycle of the noisy signal');
xlabel('Time Frame (sec)'); ylabel('Voltage (mv)');

OneCycleECGPSD=figure('Name','PSDs of Clean & Noisy ECGs from Lab 1'); % Create a new figure

CycleCleanECG = FindPSD(CycleCleanECG_Voltage, length(CycleCleanECG_Voltage), 200);
subplot(4,1,1), hold on;
plot(CycleCleanECG.HalfFs, CycleCleanECG.PSD);
plot(CycleCleanECG.HalfFs, -20, 'r'); % plot dB line
hold off;
title('PSD - 1 cycle of clean ECG (fs = 200 Hz) (N = # samples in signal)');
ylabel('Magnitude (dB)');
axis([0 (CycleCleanECG.fs/2) -50 5]); % If leave to auto it will be hard to compare

CycleNoisyECG = FindPSD(CycleNoisyECG_Voltage, length(CycleNoisyECG_Voltage), 200);
subplot(4,1,2), hold on;
plot(CycleNoisyECG.HalfFs, CycleNoisyECG.PSD);
plot(CycleNoisyECG.HalfFs, -20, 'r'); % plot dB line
hold off;
title('PSD - 1 cycle of noisy ECG (fs = 200 Hz) (N = # samples in signal)');
ylabel('Magnitude (dB)');
axis([0 (CycleNoisyECG.fs/2) -50 5]); % If leave to auto it will be hard to compare

% ---------------------------------------------------- Experiment 1 Part 2
% Obtain the PSDs of one cardiac cycle of clean and noisy ECG with 2048
% with a larger number of samples (zero padding).  Resolution in freq axis
% will increase.

CycleCleanECGZP = FindPSD(CycleCleanECG_Voltage, 2048, 200);
subplot(4,1,3), hold on;
plot(CycleCleanECGZP.HalfFs, CycleCleanECGZP.PSD);
plot(CycleCleanECGZP.HalfFs, -20, 'r'); % plot dB line
hold off;
title('PSD - 1 cycle of clean ECG (fs = 200 Hz) (N = 2048)');
ylabel('Magnitude (dB)');
axis([0 (CycleCleanECGZP.fs/2) -50 5]); % If leave to auto it will be hard to compare

CycleNoisyECGZP = FindPSD(CycleNoisyECG_Voltage, 2048, 200);
subplot(4,1,4), hold on;
plot(CycleNoisyECGZP.HalfFs, CycleNoisyECGZP.PSD);
plot(CycleNoisyECGZP.HalfFs, -20, 'r'); % plot dB line
hold off;
title('PSD - 1 cycle of noisy ECG (fs = 200 Hz) (N = 2048)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (CycleNoisyECGZP.fs/2) -50 5]); % If leave to auto it will be hard to compare

% ---------------------------------------------------- Experiment 1 Part 3
% Obtain the PSD of 20 cycles of ECG

TwentyCycleCleanECG=figure('Name','20 Cycles of Clean ECG'); % Create a new figure

TwentyCycleCleanECG_Time = CleanECG(1:2900,1);   %  copy column 1 of the matrix to Time
TwentyCycleCleanECG_Voltage = CleanECG(1:2900,2);   %  copy column 2 of the matrix to Voltage
subplot(2,1,1), plot(TwentyCycleCleanECG_Time,TwentyCycleCleanECG_Voltage);
title('ECG #1 - 20 cycles of clean signal');
xlabel('Time Frame (sec)'); ylabel('Voltage (mv)');

TwentyCycleCleanECG = FindPSD(TwentyCycleCleanECG_Voltage, length(TwentyCycleCleanECG_Voltage), 200);
subplot(2,1,2), hold on;
plot(TwentyCycleCleanECG.HalfFs, TwentyCycleCleanECG.PSD);
plot(TwentyCycleCleanECG.HalfFs, -20, 'r'); % plot db line
hold off;
title('PSD - 20 Cycles of clean ECG (fs = 200 Hz) (N = # samples in signal)');
xlabel('Frequency Hz'); ylabel('Magnitude (dB)');
axis([0 (TwentyCycleCleanECG.fs/2) -50 5]); % If leave to auto it will be hard to compare