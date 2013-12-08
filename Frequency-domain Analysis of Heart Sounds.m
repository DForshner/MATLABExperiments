% Frequency-domain Analysis of Heart Sounds

close all; clear all; clc; % Clear everything

% -------------------------------------------------------------------------
% ---------------------------------------------------------------- pec1.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec1.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec1.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec1.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec1.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% ---------------------------------------------------------------- pec1.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec1.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec1.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec1.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec1.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% --------------------------------------------------------------- pec22.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec22.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec22.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec22.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)
% Find the Sync Average of the PSDs for each segment

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec22.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% --------------------------------------------------------------- pec33.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec33.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec33.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec33.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec33.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% -------------------------------------------------------------- pec41.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec41.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec41.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec41.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec44.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% --------------------------------------------------------------- pec42.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec42.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec42.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec42.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec42.wav'); % Make an audio wave file of the PCG signal

% -------------------------------------------------------------------------
% --------------------------------------------------------------- pec52.dat
% -------------------------------------------------------------------------

% ------------------------------------------------------- Part 1: Load File
Signal = load('pec52.dat'); % loading ascii file into an array

% Separating pcg, ecg, and carotid signals
PCG = Signal(:,1);
ECG = Signal(:,2);
Carotid = Signal(:,3);
Fs = 1000; % Sampling frequency is 1000Hz
TimeAxis = [1:length(PCG)]/Fs;

% ----------------------------------------------------- Part 2: Segment PCG
DownSampledECG = DownSample(ECG, Fs); % Down sample ECG signal to 200 Hz
DownSampledTimeAxis = TimeAxis(5:5:end); % Takes every fifth sample

% Apply the Pan–Tompkins algorithmto the ECG signal and detect the QRS complexes.
QRSPeakIndex = FindQRS(DownSampledECG, DownSampledTimeAxis, Fs);

% starting from the beginning of each QRS complex.
QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs);
QRSBeginIndex = QRSBeginIndex * 5; % Rescale the sample index up to 1000 Hz

% Segment the systolic portions of the PCG signal by selecting a window of duration 300 − 350 ms
PCGS1BeginIndex = QRSBeginIndex;
PCGS1EndIndex = QRSBeginIndex + 300; % 300 ms = 350 samples

% Plot PCG with S1 segments indicated
PCGS1Segmented = figure('Name','Pec52.dat - PCG with S1 sounds segmented'); % Create a new figure
plot(TimeAxis, PCG); hold on; 
plot(TimeAxis(PCGS1BeginIndex), PCG(PCGS1BeginIndex), 'r*');
plot(TimeAxis(PCGS1EndIndex), PCG(PCGS1EndIndex), 'r+'); hold off;
title('pec52.dat - PCG with S1 sounds segmented (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('PCG (Mv)'); axis auto;

% --------------------------------------- Part 3: Sync average segment PSDs
% Find the Sync Average of the PSDs for each segment
SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)

% ----------------------------------------------- Part 4: Create Audio file
MakeWave(PCG, 'pec52.wav'); % Make an audio wave file of the PCG signal