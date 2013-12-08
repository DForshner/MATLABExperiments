% Frequency-domain Analysis of Biomedical Signals
% Experiment 3
%
% From the EMG signal you acquired in Lab 1, select one segment for each level of muscular
% contraction. Obtain and analyze the PSD of each segment.

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

load EMGAndForceVsTime.txt;  % Load file into memory
EMG_Time = EMGAndForceVsTime(:,1);   % copy column 1 of the matrix to Time
EMG_Force = EMGAndForceVsTime(:,2);   % copy column 2 of the matrix to Voltage
EMG_Voltage = EMGAndForceVsTime(:,3);   % copy column 3 of the matrix to Force

% ---------------------------------------------------------- Part 1
% Normalize the Force vector
EMG_ForceShifted = EMG_Force - min(EMG_Force);
EMG_ForceNorm = EMG_ForceShifted / max(EMG_ForceShifted) * 100; % Note that min(EMG_ForceShifted is zero)

%Subtract the mean from the EMG signal to remove any bias
EMG_Voltage = EMG_Voltage - mean(EMG_Voltage);

% ---------------------------------------------------------- Part 2
% Manually identify portions (segments) corresponding to each level of
% contraction within which the force remains almost constant
IdentSegTime = [4.0, 6.2, 9.8, 12.6, 17.1, 20.1, 24.5, 28.7, 33.9, 37.8, 42.5, 47, 52.6 ,57.1];

% Plot of the EMG (mV) and normalized force (in percent % MVC) vs  Time
% (Sec) with segments marked
EMG_VoltageForceVsTimeSeg=figure('Name','EMG signal and Normalized force vs. Time'); % Create a new figure

subplot(2,1,1); plot(EMG_Time, EMG_ForceNorm)
hold on;
plot(IdentSegTime,0,'r*');  % Highlight 
hold off;
ylabel('Normalized Force (%MVC)'); xlabel('Time (Sec)');

subplot(2,1,2); plot(EMG_Time, EMG_Voltage)
hold on;
plot(IdentSegTime,0,'r*');  % Highlight 
hold off;
ylabel('EMG Signal (mV)'); xlabel('Time (Sec)');
title('EMG signal and Normalized force vs. Time');
axis auto; 

% ---------------------------------------------------------- Part 3
% For each segment of the EMG identified obtain and analyze the PSD of each segment.

EMGFigure = figure('Name','PSDs force contraction segments of the EMG signal'); % Create a new figure

SampleInterval = 0.0005; % The EMG signals are sampled with a period of 0.0005s
Fs = 1 / SampleInterval; % Sampling freqency

for Counter=1:1:7
    Segment = Counter * 2 - 1;
    disp(['Segment = ', num2str(Segment)]); % Display in console
    disp(['Start = ', num2str(IdentSegTime(Segment))]); % Display in console
    disp(['Stop = ', num2str(IdentSegTime(Segment + 1))]); % Display in console
    
    % Identify the segments of interest from the EMG signals
    s1 = IdentSegTime(Segment) / SampleInterval; % First point index
    s2 = IdentSegTime(Segment + 1) / SampleInterval; % Second point index
    EMG_ForceSeg = EMG_ForceNorm(s1:s2); % Segment of interest of force signal
    EMG_VoltageSeg = EMG_Voltage(s1:s2); % Segment of interest of voltage signal  
        
    EMG = FindPSD(EMG_VoltageSeg, length(EMG_VoltageSeg), (Fs));
    subplot(7,1,Counter), hold on;
    plot(EMG.HalfFs, EMG.PSD);
    plot(EMG.HalfFs, -15, 'r'); % plot db line
    hold off;
    title(['PSD of seg #' num2str(Counter),' (Mean Force = ', num2str(mean(EMG_ForceSeg)),' %MCV) (Fs = 2000 Hz) (N = # samples)']);
    ylabel('Magnitude (dB)');
    axis([0 (EMG.fs/2) -40 5]); % If leave to auto it will be hard to compare   
end

xlabel('Frequency (Hz)');