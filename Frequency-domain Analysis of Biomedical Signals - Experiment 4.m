% Frequency-domain Analysis of Biomedical Signals
% Experiment 1
%
% Part 1 - Using the program pcg3read.m, read the dataﬁle pec1.dat; the ﬁles 
% are available at The program gives a plot of three channels: ECG, PCG, 
% and carotid pulse. The sampling rate is fs = 1000 Hz per channel.
%
% Part 2 - Select a systolic segment of the PCG signal and obtain its
% PSD. Select a diastolic segment and obtain its PSD.
%
% Part 3 - Repeat the above with the data in the ﬁle pec33.dat. This signal 
% contains systolic systolic and diastolic segments. Subtract 52 ms from the 
% dicrotic notch position to obtain the initial time instant of the diastolic 
% segment; Obtain the PSDs of the systolic and diastolic segments.

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

% ---------------------------------------------------- Experiment 4 Part 1
% Read and plot the combined PCG/ECG/carotid pulse datafile

load 'pec1.dat'; % Load the PCG, ECG, and carotid pulse data

fs = 1000; % Sampling frequency is fs = 1000 Hz
PCGSignal = pec1(:,1); % PCG
ECGSignal = pec1(:,2); % ECG
CarotidSignal = pec1(:,3); % carotid 
clear pec1;

N = length(PCGSignal); % number of samples
T = 1/fs; % period
t = [0 : (N - 1)] * T; % time axis

%Plot the signals
Pec1Plot = figure('Name','Pec1.dat - PCG, ECG, and carotid pulse signals'); % Create a new figure
subplot(3,1,1); plot(t, PCGSignal);
title('Pec1.dat - PCG, ECG, and carotid pulse signals');
ylabel('PCG (Unknown)');
axis tight;
subplot(3,1,2); plot(t, ECGSignal);
ylabel('ECG (Unknown)');
axis tight;
subplot(3,1,3); plot(t, CarotidSignal);
ylabel('Carotid pulse (Unknown)'); xlabel('Time (Sec)');
axis tight;

% ---------------------------------------------------- Experiment 4 Part 2
% Select a systolic segment of the PCG signal and obtain its
% PSD. Select a diastolic segment and obtain its PSD.

% Manually identify the systolic segment & diastolic segment of the PCG
IdentSegSample = [4089, 4455, 5053];

subplot(3,1,1); hold on;
plot((IdentSegSample*T),PCGSignal(IdentSegSample),'r*');  % Highlight 
hold off;
xlim([4 5.5]);
subplot(3,1,2); hold on;
plot((IdentSegSample*T),ECGSignal(IdentSegSample),'r*');  % Highlight 
hold off;
xlim([4 5.5]);
subplot(3,1,3); hold on;
plot((IdentSegSample*T),CarotidSignal(IdentSegSample),'r*');  % Highlight 
hold off;
xlim([4 5.5]);

% Find & plot the PSD of the systolic PCG signal
Pec1PCGPSDPlot = figure('Name','Pec1.dat - PSDs for systolic & diastolic PCG signal'); % Create a new figure

PCGSysSeg = PCGSignal(IdentSegSample(1):IdentSegSample(2)); % Segment of interest
subplot(4,1,1), plot( t(IdentSegSample(1):IdentSegSample(2)) , PCGSysSeg);
title('Pec1.dat - Systolic segment of PCG signal');
ylabel('PCG (Unknown)'); xlabel('Time (Sec)');
axis auto;

PCGSystolic = FindPSD(PCGSysSeg, 2048, fs);
subplot(4,1,3), hold on;
plot(PCGSystolic.HalfFs, PCGSystolic.PSD);
plot(PCGSystolic.HalfFs, -20, 'r');  % plot the -20db line
hold off;
title(['PSD of Systolic PCG segment of PCG (Fs = 1000 Hz) (N = # 2048)']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (PCGSystolic.fs/2) -50 5]); % If leave to auto it will be hard to compare  

% Find the PSD of the diastolic PCG signal
PCGDiaSeg = PCGSignal(IdentSegSample(2):IdentSegSample(3)); % Segment of interest

subplot(4,1,2), plot( t(IdentSegSample(2):IdentSegSample(3)) , PCGDiaSeg);

title(['Pec1.dat - Diastolic segment of PCG signal']);
ylabel('PCG (Unknown)'); xlabel('Time (Sec)');
axis auto;

PCGDiastolic = FindPSD(PCGDiaSeg, 2048, fs);
subplot(4,1,4), hold on; 
plot(PCGDiastolic.HalfFs, PCGDiastolic.PSD);
plot(PCGDiastolic.HalfFs, -20, 'r');  % plot the -20db line
hold off;
title(['PSD of Diastolic PCG segment of PCG (Fs = 1000 Hz) (N = # 2048)']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (PCGDiastolic.fs/2) -50 5]); % If leave to auto it will be hard to compare

% ---------------------------------------------------- Experiment 4 Part 3
% Repeat the above with the data in the ﬁle pec33.dat. Subtract 52 ms from the 
% dicrotic notch position to obtain the initial time instant of the diastolic 
% segment; Obtain the PSDs of the systolic and diastolic segments.

load 'pec33.dat'; % Load the PCG, ECG, and carotid pulse data

fs = 1000; % Sampling frequency is fs = 1000 Hz
PCGSignal = pec33(:,1); % PCG
ECGSignal = pec33(:,2); % ECG
CarotidSignal = pec33(:,3); % carotid 
clear pec33;

N = length(PCGSignal); % number of samples
T = 1/fs; % period
t = [0 : (N - 1)] * T; % time axis

%Plot the signals
Pec33Plot = figure('Name','Pec33.dat - PCG, ECG, and carotid pulse signals'); % Create a new figure
subplot(3,1,1); plot(t, PCGSignal);
title('PCG, ECG, and carotid pulse signals');
ylabel('PCG (Unknown)');
axis tight;
subplot(3,1,2); plot(t, ECGSignal);
ylabel('ECG (Unknown)');
axis tight;
subplot(3,1,3); plot(t, CarotidSignal);
ylabel('Carotid pulse (Unknown)'); xlabel('Time (Sec)');
axis tight;

% Manually identify the systolic segment & diastolic segment of the PCG
IdentSegSample = [4125, 4449, 4683];

subplot(3,1,1); hold on;
plot((IdentSegSample(1)*T),PCGSignal(IdentSegSample(1)),'r*');  % Highlight
plot((IdentSegSample(3)*T),PCGSignal(IdentSegSample(3)),'r*');  % Highlight
hold off;
xlim([4 5]);
subplot(3,1,2); hold on;
plot((IdentSegSample*T),ECGSignal(IdentSegSample),'r*');  % Highlight 
hold off;
xlim([4 5]);
subplot(3,1,3); hold on;
plot((IdentSegSample*T),CarotidSignal(IdentSegSample),'r*');  % Highlight 
hold off;
xlim([4 5]);

% Subtract 52 ms from the dicrotic notch position to obtain the initial time 
% instant of the diastolic segment
IndexShift = round(0.052 / T); % determine number of samples to shift 0.052s
IdentSegSample(2) = IdentSegSample(2) - IndexShift; % shift the diastolic seg

subplot(3,1,1); hold on;
plot((IdentSegSample(2)*T),PCGSignal(IdentSegSample(2)),'r+');  % Highlight 
hold off;

% Find & plot the PSD of the systolic PCG signal
Pec33PCGPSDPlot = figure('Name','Pec33.dat - PSDs for systolic & diastolic PCG signal'); % Create a new figure

PCGSysSeg = PCGSignal(IdentSegSample(1):IdentSegSample(2)); % Segment of interest
subplot(4,1,1), plot( t(IdentSegSample(1):IdentSegSample(2)) , PCGSysSeg);
title('Pec33.dat - Systolic segment of PCG signal');
ylabel('PCG (Unknown)'); xlabel('Time (Sec)');
axis auto;

PCGSystolic = FindPSD(PCGSysSeg, 2048, fs);
subplot(4,1,3), hold on;
plot(PCGSystolic.HalfFs, PCGSystolic.PSD);
plot(PCGSystolic.HalfFs, -20, 'r');  % plot the -20db line
hold off;
title('PSD of Systolic PCG segment of PCG (Fs = 1000 Hz) (N = # 2048)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (PCGSystolic.fs/2) -50 5]); % If leave to auto it will be hard to compare  

% Find the PSD of the diastolic PCG signal
PCGDiaSeg = PCGSignal(IdentSegSample(2):IdentSegSample(3)); % Segment of interest

subplot(4,1,2), plot( t(IdentSegSample(2):IdentSegSample(3)) , PCGDiaSeg);
title('Pec33.dat - Diastolic segment of PCG signal');
ylabel('PCG (Unknown)'); xlabel('Time (Sec)');
axis auto;

PCGDiastolic = FindPSD(PCGDiaSeg, 2048, fs);
subplot(4,1,4), hold on;
plot(PCGDiastolic.HalfFs, PCGDiastolic.PSD);
plot(PCGDiastolic.HalfFs, -20, 'r'); % plot the -20db line
hold off;
title('PSD of Diastolic PCG segment of PCG (Fs = 1000 Hz) (N = # 2048)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([0 (PCGDiastolic.fs/2) -50 5]); % If leave to auto it will be hard to compare
